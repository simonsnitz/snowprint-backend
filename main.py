from src.blast import blast
from src.get_genome_coordinates import get_genome_coordinates, get_genome_coordinates_batch
from src.accID2operon import acc2operon
from src.fetch_promoter import fetch_promoter
from src.fetch_operator import fetch_operator
from src.troubleshoot import troubleshoot

import re
import pandas as pd
import requests
import json
import sys
import ast
import os
import base64
import boto3
import time
import hashlib

from decimal import Decimal

def set_params(param_obj):
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

    return blast_params, promoter_params, operator_params

# (1) BLAST protein. return a dataframe

def perform_blast(blast_params, promoter_params, operator_params, data):
    genome_choice = data["genomeChoice"]
    filter_redundant = data["filter"]
    acc = data["acc"]
    input_method = data["method"]
    max_homologs = data["homologs"]

    blast_df = blast(acc, input_method, blast_params, max_seqs=500)
    homolog_dict = []


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
            homolog_dict = filter_blastDf(blast_df)
        else:
            homolog_dict = [
                {"Uniprot Id": row["Uniprot Id"], "identity": row["Identity"],"coverage": row["Coverage"]}
                for i, row in blast_df.iterrows()
                ]

        # limit search to specified number of homologs
        homolog_dict = homolog_dict[0:max_homologs]
        
        ### DISPLAY homolog_dict in the frontend
        print("blast finished.")

    # (2) Get genome coordianates. Return a dataframe

    if genome_choice == "batch":
        homolog_dict = get_genome_coordinates_batch(homolog_dict)

        #TODO: I get an error here sometimes.
        if homolog_dict == None:
            print("Failed fetching genome coordinates. Try fetching these individually (advanced options)")
        else:
            homolog_dict = [i for i in homolog_dict if i != None]


    elif genome_choice == "individually":

        updated_homolog_dict = []
        for i in range(0, len(homolog_dict)):
            updated_homolog_dict.append(get_genome_coordinates(homolog_dict[i]))

        # Remove entries without any genome coordinates
        homolog_dict = [i for i in updated_homolog_dict if i != None]


    homolog_dict = [i for i in homolog_dict if "Genome" in i.keys()]
    coordinates_df = pd.DataFrame(homolog_dict).drop(columns=["identity", "coverage"])

    ### DISPLAY coordinates_df in the frontend
    print("genome coordinates fetched.")


    # (3) Extract predicted operators for each homolog. return a dataframe

    for i in range(0, len(homolog_dict)):
        homolog_dict[i]["operon"] = acc2operon(homolog_dict[i])
        # Deal with cases where operon fetching fails
        try:
            homolog_dict[i]["promoter"] = fetch_promoter(homolog_dict[i]["operon"], promoter_params)
        except:
            homolog_dict[i]["promoter"] = None


    operator_dict = fetch_operator(homolog_dict, operator_params)
    operator_df = pd.DataFrame(operator_dict["aligned_seqs"])

    ### DISPLAY operator dataframe in the frontend.
    print("operators fetched.")


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

    primary = passed_in_data["acc"]
    coordinates = return_data["coordinates_df"].to_dict('records')
    aligned_seqs = return_data["aligned_seqs"].to_dict('records')

    # Construct a sha-256 hash for future reference if these advanced options have been done before
    json_string = json.dumps(passed_in_data, sort_keys=True)
    hashed = hashlib.sha256(json_string.encode()).hexdigest()

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

if __name__ == "__main__":

    # The following code may look strange but there is a method to the madness:
    # My goal is to create ECS tasks that are completely async - they get data and then write to a DynamoDB when they are done
    # I didn't want them to be dependent on any other AWS service to further complicate things
    # This task is started by lambda to keep costs as low as possible
    # The challenge became - how do we get JSON data into an ECS container?
    # One solution could be using S3 but then we'd need to create unique JSON files in S3 and feed that specific file name to the container
    # ^ This felt like a waste of a service to just hold a json object
    # Another idea was just to create a file in /tmp/ in lambda but there's no way to specify a volume source location in boto3
    # ^ You can configure a mountPoint in the ECS task definition but again, using NamedTemporaryFile results in a unique name which can't be overwritten in containerOverrides
    # Then - I tried to just stringify the JSON object and pass it into the container with -e through docker
    # This would have worked but an object in the format of '{"some":"value"}' threw an error since docker through it was the repository name no matter what way to shook it
    
    # Thus - this solution is born:
    # 1. Lambda recieves the event body from the API requests
    # 2. It encodes it into base64 (this format ensures docker does not think it's a repository name) and passes in the environment variable (no need for a file)
    # 3. This container reads the environment variable
    # 4. It then must decode the base64, decode the bytestring, then double json.loads to get a valid python dict

    # Get environment
    PASSED_IN_JSON = os.environ['PASSED_IN_JSON']
    TABLE_NAME = os.environ['TABLE_NAME']
    PASSED_UUID = os.environ['UUID']

    # Decode base64 JSON
    decode_base64 = base64.b64decode(PASSED_IN_JSON)
    decode_bytes = decode_base64.decode('utf-8')
    data = json.loads(decode_bytes)

    # Start snowprint
    blast_params, promoter_params, operator_params = set_params(data)
    return_data = perform_blast(blast_params, promoter_params, operator_params, data)

    # Write results
    extracted_dict = format_homologs(return_data["homolog_dict"])
    write_results_to_db(TABLE_NAME, extracted_dict, return_data, data, PASSED_UUID)
