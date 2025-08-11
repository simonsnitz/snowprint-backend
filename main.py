from src.blast import blast
from src.get_genome_coordinates import get_genome_coordinates, get_genome_coordinates_batch
from src.accID2operon import acc2operon
from src.fetch_promoter import fetch_promoter
from src.fetch_operator import fetch_operator
from src.troubleshoot import troubleshoot
from src.circuit_breaker import embl_api_call, uniprot_api_call, blast_call

import re
import pandas as pd
import requests
import json
import sys
import signal
import ast
import os
import base64
import boto3
import time
import hashlib
import logging
import psutil
import threading
from datetime import datetime, timezone

from decimal import Decimal

# Configure structured logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)
logger = logging.getLogger(__name__)

# Global variables for resource monitoring
start_time = None
request_id = None
task_uuid = None
resource_monitor_thread = None
stop_monitoring = threading.Event()
timeout_thread = None

# Maximum execution time (45 minutes)
MAX_EXECUTION_TIME = 45 * 60

# Resource monitoring function
def monitor_resources():
    """Monitor CPU and memory usage during execution"""
    process = psutil.Process()
    max_memory = 0
    max_cpu = 0
    
    while not stop_monitoring.is_set():
        try:
            memory_mb = process.memory_info().rss / 1024 / 1024
            cpu_percent = process.cpu_percent()
            
            max_memory = max(max_memory, memory_mb)
            max_cpu = max(max_cpu, cpu_percent)
            
            # Log every 30 seconds
            if int(time.time()) % 30 == 0:
                logger.info(f"[{request_id}] Resource usage", extra={
                    'task_uuid': task_uuid,
                    'memory_mb': memory_mb,
                    'cpu_percent': cpu_percent,
                    'max_memory_mb': max_memory,
                    'max_cpu_percent': max_cpu,
                    'runtime_seconds': time.time() - start_time if start_time else 0
                })
            
            time.sleep(1)
        except Exception as e:
            logger.error(f"[{request_id}] Error monitoring resources: {str(e)}")
            break

def timeout_handler():
    """Monitor total execution time and force exit if exceeded"""
    global start_time
    
    while not stop_monitoring.is_set():
        if start_time and (time.time() - start_time) > MAX_EXECUTION_TIME:
            logger.critical(f"[{request_id}] Maximum execution time ({MAX_EXECUTION_TIME}s) exceeded - forcing shutdown", extra={
                'task_uuid': task_uuid,
                'runtime_seconds': time.time() - start_time,
                'max_allowed_seconds': MAX_EXECUTION_TIME
            })
            
            # Write timeout error to database
            try:
                central_error_handling(f"Task timeout after {MAX_EXECUTION_TIME} seconds", {})
            except Exception:
                pass
            
            # Force exit
            os._exit(1)
            
        time.sleep(30)  # Check every 30 seconds

# Graceful shutdown function
def signal_handler(signal_num, frame):
    logger.warning(f"[{request_id}] Graceful shutdown activated", extra={
        'task_uuid': task_uuid,
        'signal': signal_num,
        'runtime_seconds': time.time() - start_time if start_time else 0
    })
    
    # Stop resource monitoring
    stop_monitoring.set()
    if resource_monitor_thread and resource_monitor_thread.is_alive():
        resource_monitor_thread.join(timeout=5)
    
    # Log final resource usage
    try:
        process = psutil.Process()
        memory_mb = process.memory_info().rss / 1024 / 1024
        logger.info(f"[{request_id}] Final resource usage", extra={
            'task_uuid': task_uuid,
            'final_memory_mb': memory_mb,
            'total_runtime_seconds': time.time() - start_time if start_time else 0
        })
    except Exception:
        pass

# Map signals to handler
signal.signal(signal.SIGTERM, signal_handler)
signal.signal(signal.SIGINT, signal_handler)

def set_params(param_obj):
    try:
        blast_params = {
            "ident_cutoff": param_obj["identity"],
            "cov_cutoff": param_obj["coverage"]
        }

        promoter_params = {
            "min_length": param_obj["minLength"],
            "max_length": param_obj["maxLength"]
        }

        operator_params = {
            "extension_length": param_obj["extension"],
            "win_score": param_obj["alignMatch"],
            "loss_score": param_obj["alignMismatch"],
            "spacer_penalty": param_obj["penalty"],
            "gap_open": param_obj["gapOpen"],
            "gap_extend": param_obj["gapExtend"],
            "align_match": param_obj["alignMatch"],
            "align_mismatch": param_obj["alignMismatch"],
            "min_operator_length": param_obj["minOperator"],
            "max_operator_length": param_obj["maxOperator"],
            "seq_to_align": param_obj["seqToAlign"],
            "search_method": param_obj["conservation"]
        }
    except Exception as e:
        print(e)
        raise Exception('Set params failed')

    return blast_params, promoter_params, operator_params

# (1) BLAST protein. return a dataframe

def perform_blast(blast_params, promoter_params, operator_params, data):
    stage_start = time.time()
    genome_choice = data["genomeChoice"]
    filter_redundant = data["filter"]
    acc = data["acc"]
    input_method = data["method"]
    max_homologs = data["homologs"]
    homolog_dict = []

    logger.info(f"[{request_id}] Starting BLAST stage", extra={
        'task_uuid': task_uuid,
        'stage': 'blast',
        'input_method': input_method,
        'acc': acc,
        'max_homologs': max_homologs,
        'params': blast_params
    })

    try:
        # Use circuit breaker for BLAST call
        blast_df = blast_call(blast, acc, input_method, blast_params, max_seqs=500)
        blast_duration = time.time() - stage_start
        
        logger.info(f"[{request_id}] BLAST completed", extra={
            'task_uuid': task_uuid,
            'stage': 'blast',
            'duration_seconds': blast_duration,
            'results_count': len(blast_df) if not blast_df.empty else 0
        })
        
    except Exception as e:
        blast_duration = time.time() - stage_start
        logger.error(f"[{request_id}] BLAST failed", extra={
            'task_uuid': task_uuid,
            'stage': 'blast',
            'duration_seconds': blast_duration,
            'error': str(e)
        })
        raise Exception(e)

    if not blast_df.empty:

        #if 'filter redundant' box checked, filter out homologs that have the same %identity and %coverage
        def filter_blastDf(blast_df):
            
            ident_covs = []
            for i, row in blast_df.iterrows():
                entry =   {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
                to_compare =   {"identity": row["Identity"],"coverage": row["Coverage"]}
                if to_compare not in ident_covs:
                    homolog_dict.append(entry)
                    ident_covs.append(to_compare)
            return homolog_dict

        if filter_redundant:
            try:
                homolog_dict = filter_blastDf(blast_df)
            except Exception as e:
                print('Error trying to filter blastDf')
                raise Exception(e)
        else:
            homolog_dict = [
                {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
                for i, row in blast_df.iterrows()
                ]

        # limit search to specified number of homologs
        homolog_dict = homolog_dict[0:max_homologs]
        
        logger.info(f"[{request_id}] BLAST homologs processed", extra={
            'task_uuid': task_uuid,
            'stage': 'blast_filtering',
            'filtered_homologs': len(homolog_dict),
            'filter_redundant': filter_redundant
        })

    # (2) Get genome coordinates. Return a dataframe
    coordinates_start = time.time()
    
    logger.info(f"[{request_id}] Starting genome coordinates stage", extra={
        'task_uuid': task_uuid,
        'stage': 'coordinates',
        'genome_choice': genome_choice,
        'homologs_to_process': len(homolog_dict)
    })

    if genome_choice == "batch":
        try:
            # Use circuit breaker for batch genome coordinates
            homolog_dict = embl_api_call(get_genome_coordinates_batch, homolog_dict)
            coordinates_duration = time.time() - coordinates_start
            
            logger.info(f"[{request_id}] Batch genome coordinates completed", extra={
                'task_uuid': task_uuid,
                'stage': 'coordinates',
                'duration_seconds': coordinates_duration,
                'successful_coordinates': len([h for h in homolog_dict if h and "Genome" in h.keys()]) if homolog_dict else 0
            })
            
        except Exception as e:
            coordinates_duration = time.time() - coordinates_start
            logger.error(f"[{request_id}] Batch genome coordinates failed", extra={
                'task_uuid': task_uuid,
                'stage': 'coordinates',
                'duration_seconds': coordinates_duration,
                'error': str(e)
            })
            raise Exception(e)

        if homolog_dict == None:
            error_msg = "Failed fetching genome coordinates. Try fetching these individually (advanced options)"
            logger.error(f"[{request_id}] {error_msg}", extra={'task_uuid': task_uuid, 'stage': 'coordinates'})
            raise Exception(error_msg)
        else:
            homolog_dict = [i for i in homolog_dict if i != None]

    elif genome_choice == "individually":
        updated_homolog_dict = []
        api_call_count = 0
        failed_calls = 0
        
        for i in range(0, len(homolog_dict)):
            try:
                api_call_count += 1
                # Use circuit breaker for individual genome coordinates
                result = embl_api_call(get_genome_coordinates, homolog_dict[i])
                updated_homolog_dict.append(result)
                
                if (i + 1) % 5 == 0:  # Log every 5 calls
                    logger.info(f"[{request_id}] Individual coordinates progress", extra={
                        'task_uuid': task_uuid,
                        'stage': 'coordinates_individual',
                        'processed': i + 1,
                        'total': len(homolog_dict),
                        'api_calls': api_call_count,
                        'failed_calls': failed_calls
                    })
                    
            except Exception as e:
                failed_calls += 1
                logger.warning(f"[{request_id}] Individual coordinate fetch failed", extra={
                    'task_uuid': task_uuid,
                    'stage': 'coordinates_individual',
                    'homolog_index': i,
                    'error': str(e)
                })
                raise Exception(e)

        coordinates_duration = time.time() - coordinates_start
        logger.info(f"[{request_id}] Individual genome coordinates completed", extra={
            'task_uuid': task_uuid,
            'stage': 'coordinates',
            'duration_seconds': coordinates_duration,
            'total_api_calls': api_call_count,
            'failed_calls': failed_calls
        })

        # Remove entries without any genome coordinates
        homolog_dict = [i for i in updated_homolog_dict if i != None]

    homolog_dict = [i for i in homolog_dict if "Genome" in i.keys()]
    coordinates_df = pd.DataFrame(homolog_dict).drop(columns=["identity", "coverage"])

    logger.info(f"[{request_id}] Genome coordinates stage completed", extra={
        'task_uuid': task_uuid,
        'stage': 'coordinates',
        'final_homologs_with_genome': len(homolog_dict)
    })


    # (3) Extract predicted operators for each homolog. return a dataframe
    operator_start = time.time()
    
    logger.info(f"[{request_id}] Starting operator extraction stage", extra={
        'task_uuid': task_uuid,
        'stage': 'operators',
        'homologs_to_process': len(homolog_dict)
    })

    promoter_fetch_count = 0
    promoter_failures = 0
    
    for i in range(0, len(homolog_dict)):
        homolog_dict[i]["operon"] = acc2operon(homolog_dict[i])
        # Deal with cases where operon fetching fails
        try:
            homolog_dict[i]["promoter"] = fetch_promoter(homolog_dict[i]["operon"], promoter_params)
            promoter_fetch_count += 1
        except Exception as e:
            homolog_dict[i]["promoter"] = None
            promoter_failures += 1
            logger.debug(f"[{request_id}] Promoter fetch failed for homolog {i}: {str(e)}")

    logger.info(f"[{request_id}] Promoter fetching completed", extra={
        'task_uuid': task_uuid,
        'stage': 'promoters',
        'successful_promoters': promoter_fetch_count,
        'failed_promoters': promoter_failures,
        'total_homologs': len(homolog_dict)
    })

    operator_dict = fetch_operator(homolog_dict, operator_params)
    operator_df = pd.DataFrame(operator_dict["aligned_seqs"])
    
    operator_duration = time.time() - operator_start
    
    logger.info(f"[{request_id}] Operator extraction completed", extra={
        'task_uuid': task_uuid,
        'stage': 'operators',
        'duration_seconds': operator_duration,
        'aligned_sequences': len(operator_dict["aligned_seqs"]) if operator_dict["aligned_seqs"] else 0,
        'consensus_score': operator_dict.get("consensus_score", "unknown"),
        'num_seqs': operator_dict.get("num_seqs", "unknown")
    })


    # (4) Process results

    # Display metrics
    metric1 = operator_dict["consensus_score"]
    metric2 = operator_dict["num_seqs"]

    # Show where the predicted operator is located within the native promoter
    for i in homolog_dict:
        if i["promoter"]:
            [before, after] = re.split(re.escape((operator_dict["native_operator"]).upper()), i["promoter"])
            html = "<span style='color: black;'>"+str(before)+"</span>"
            html += "<span style='color: red; font-size: 16px'>"+str(operator_dict["native_operator"])+"</span>"
            html += "<span style='color: black;'>"+str(after)+"</span>"
            # DISPLAY this html code
            break

    # Display the consensus sequence
    consensus_seq = operator_dict["consensus_seq"]

    # Create & Display the consensus motif logo
    motif = operator_dict["motif"]
    motif_html = []
    color_key = {"A":"red", "T": "green", "C": "blue", "G": "#fcba03"}
    for i in motif:
        motif_html.append(
            {
                "color": str(color_key[i["base"].upper()]),
                "fontSize": 400,
                "fontWeight": 550,
                "display": 'inline-block',
                "translateY": str(1.25-i["score"]**3),
                "scaleY": str(3*i["score"]**3),
                "base": str(i["base"])
            }
        )

    return_data = {
        "homolog_dict": homolog_dict,
        "coordinates_df": coordinates_df,
        "aligned_seqs": pd.DataFrame(operator_dict["aligned_seqs"]),
        "consensus_score": operator_dict["consensus_score"],
        "num_seqs": operator_dict["num_seqs"],
        "html": html,
        "consensus_seq": consensus_seq,
        "motif_html": motif_html
    }

    return return_data
    
def format_homologs(homolog_dict):
    extracted_dict = [{"coverage": d["coverage"], "identity": d["identity"], "Uniprot Id": d["Uniprot Id"], "promoter": d["promoter"]} for d in homolog_dict]
    return extracted_dict

# def write_results_to_db(table, extracted_dict, coordinates, aligned, consensus, num, passed_in_data, PASSED_UUID):
def write_results_to_db(table, extracted_dict, return_data, passed_in_data, PASSED_UUID):
    coordinates = return_data["coordinates_df"].to_dict('records')
    aligned_seqs = return_data["aligned_seqs"].to_dict('records')

    Item={
        'PK': primary,
        'SK': PASSED_UUID,
        'homolog': extracted_dict,
        'coordinates': coordinates,
        'aligned_seq': aligned_seqs,
        'consensus_score': return_data["consensus_score"],
        'num_seqs': return_data["num_seqs"],
        "html": return_data["html"],
        "consensus_seq": return_data["consensus_seq"],
        "motif_html": return_data["motif_html"],
        'hash': hashed,
        'status': 'complete'
        }

    # Prepare DynamoDB
    parsed_data = json.loads(json.dumps(Item), parse_float=Decimal)
    dynamodb = boto3.resource('dynamodb', region_name='us-east-2')
    table = dynamodb.Table(TABLE_NAME)

    # Write data
    try:
        table.put_item(
            Item=parsed_data
        )
    except Exception as e:
        print(e)

# This function propagates error up throughout the whole application and writes it to the DB
def central_error_handling(reason, passed_in_data):
    print(f'Snowprint failure: {reason}')

    Item={
        'PK': primary,
        'SK': PASSED_UUID,
        'hash': hashed,
        'status': 'error',
        'reason': str(reason)
        }

    # Prepare DynamoDB
    parsed_data = json.loads(json.dumps(Item))
    dynamodb = boto3.resource('dynamodb', region_name='us-east-2')
    table = dynamodb.Table(TABLE_NAME)

    # Write data
    try:
        table.put_item(
            Item=parsed_data
        )
    except Exception as e:
        print(e)

    exit()

if __name__ == "__main__":
    global hashed, primary
    
    start_time = time.time()
    
    # Get environment variables
    PASSED_IN_JSON = os.environ['PASSED_IN_JSON']
    TABLE_NAME = os.environ['TABLE_NAME']
    PASSED_UUID = os.environ['UUID']
    request_id = os.environ.get('REQUEST_ID', 'unknown')
    
    task_uuid = PASSED_UUID

    # Decode base64 JSON
    decode_base64 = base64.b64decode(PASSED_IN_JSON)
    decode_bytes = decode_base64.decode('utf-8')
    data = json.loads(decode_bytes)

    # Set globals
    primary = data["acc"]
    json_string = json.dumps(data, sort_keys=True)
    hashed = hashlib.sha256(json_string.encode()).hexdigest()
    
    logger.info(f"[{request_id}] ECS task started", extra={
        'task_uuid': task_uuid,
        'primary_key': primary,
        'method': data.get('method', 'unknown'),
        'genome_choice': data.get('genomeChoice', 'unknown'),
        'homologs': data.get('homologs', 'unknown'),
        'start_time': datetime.now(timezone.utc).isoformat()
    })
    
    # Start resource monitoring
    resource_monitor_thread = threading.Thread(target=monitor_resources, daemon=True)
    resource_monitor_thread.start()
    
    # Start timeout monitoring
    timeout_thread = threading.Thread(target=timeout_handler, daemon=True)
    timeout_thread.start()

    # Start snowprint
    try:
        blast_params, promoter_params, operator_params = set_params(data)
        logger.info(f"[{request_id}] Parameters set successfully", extra={
            'task_uuid': task_uuid,
            'stage': 'setup'
        })
    except Exception as e:
        logger.error(f"[{request_id}] Parameter setup failed", extra={
            'task_uuid': task_uuid,
            'stage': 'setup',
            'error': str(e)
        })
        central_error_handling(e, data)
    
    try:
        return_data = perform_blast(blast_params, promoter_params, operator_params, data)
        
        total_runtime = time.time() - start_time
        logger.info(f"[{request_id}] All processing stages completed", extra={
            'task_uuid': task_uuid,
            'stage': 'completion',
            'total_runtime_seconds': total_runtime,
            'homologs_processed': len(return_data["homolog_dict"]),
            'consensus_score': return_data.get("consensus_score", "unknown")
        })
        
    except Exception as e:
        total_runtime = time.time() - start_time
        logger.error(f"[{request_id}] Processing failed", extra={
            'task_uuid': task_uuid,
            'stage': 'processing',
            'total_runtime_seconds': total_runtime,
            'error': str(e)
        })
        central_error_handling(e, data)

    # Write results to database
    try:
        extracted_dict = format_homologs(return_data["homolog_dict"])
        write_results_to_db(TABLE_NAME, extracted_dict, return_data, data, PASSED_UUID)
        
        total_runtime = time.time() - start_time
        logger.info(f"[{request_id}] Task completed successfully", extra={
            'task_uuid': task_uuid,
            'stage': 'final',
            'total_runtime_seconds': total_runtime,
            'final_homologs': len(extracted_dict),
            'status': 'complete'
        })
        
    except Exception as e:
        logger.error(f"[{request_id}] Database write failed", extra={
            'task_uuid': task_uuid,
            'stage': 'database',
            'error': str(e)
        })
        central_error_handling(e, data)
    
    # Stop monitoring
    stop_monitoring.set()
    if resource_monitor_thread and resource_monitor_thread.is_alive():
        resource_monitor_thread.join(timeout=5)
