"""
Circuit breaker and timeout utilities for external API calls
"""
import time
import requests
import logging
from functools import wraps
from typing import Any, Callable, Optional, Dict
from enum import Enum

logger = logging.getLogger(__name__)

class CircuitState(Enum):
    CLOSED = "CLOSED"
    OPEN = "OPEN" 
    HALF_OPEN = "HALF_OPEN"

class CircuitBreaker:
    """Simple circuit breaker implementation for API calls"""
    
    def __init__(self, failure_threshold: int = 5, reset_timeout: int = 60):
        self.failure_threshold = failure_threshold
        self.reset_timeout = reset_timeout
        self.failure_count = 0
        self.last_failure_time = None
        self.state = CircuitState.CLOSED
        self.success_count = 0
        
    def _should_attempt_reset(self) -> bool:
        """Check if we should attempt to reset from OPEN to HALF_OPEN"""
        return (
            self.state == CircuitState.OPEN and 
            self.last_failure_time and
            time.time() - self.last_failure_time > self.reset_timeout
        )
    
    def call(self, func: Callable, *args, **kwargs) -> Any:
        """Execute function with circuit breaker protection"""
        
        if self.state == CircuitState.OPEN:
            if self._should_attempt_reset():
                self.state = CircuitState.HALF_OPEN
                logger.info("Circuit breaker: Attempting reset to HALF_OPEN")
            else:
                raise Exception(f"Circuit breaker OPEN - blocking call to {func.__name__}")
        
        try:
            result = func(*args, **kwargs)
            
            # Success - reset failure count and ensure circuit is closed
            if self.state == CircuitState.HALF_OPEN:
                logger.info("Circuit breaker: Reset to CLOSED after successful call")
                self.state = CircuitState.CLOSED
                self.failure_count = 0
                
            return result
            
        except Exception as e:
            self._record_failure()
            logger.warning(f"Circuit breaker: Call to {func.__name__} failed - {str(e)}")
            raise
    
    def _record_failure(self):
        """Record a failure and potentially open the circuit"""
        self.failure_count += 1
        self.last_failure_time = time.time()
        
        if self.failure_count >= self.failure_threshold:
            self.state = CircuitState.OPEN
            logger.error(f"Circuit breaker: OPENED after {self.failure_count} failures")

# Global circuit breakers for different API endpoints
_circuit_breakers: Dict[str, CircuitBreaker] = {}

def get_circuit_breaker(name: str, failure_threshold: int = 3, reset_timeout: int = 60) -> CircuitBreaker:
    """Get or create a circuit breaker for a named endpoint"""
    if name not in _circuit_breakers:
        _circuit_breakers[name] = CircuitBreaker(failure_threshold, reset_timeout)
    return _circuit_breakers[name]

def with_timeout_and_circuit_breaker(
    timeout: int = 30,
    circuit_breaker_name: Optional[str] = None,
    failure_threshold: int = 3,
    reset_timeout: int = 60
):
    """Decorator to add timeout and circuit breaker protection to functions"""
    
    def decorator(func: Callable) -> Callable:
        @wraps(func)
        def wrapper(*args, **kwargs):
            # Set up timeout for requests
            if 'timeout' not in kwargs and hasattr(requests, 'get'):
                kwargs['timeout'] = timeout
            
            # Use circuit breaker if specified
            if circuit_breaker_name:
                breaker = get_circuit_breaker(circuit_breaker_name, failure_threshold, reset_timeout)
                return breaker.call(func, *args, **kwargs)
            else:
                return func(*args, **kwargs)
                
        return wrapper
    return decorator

def safe_api_call(
    func: Callable,
    *args,
    timeout: int = 30,
    max_retries: int = 2,
    circuit_breaker_name: Optional[str] = None,
    **kwargs
) -> Any:
    """
    Make a safe API call with timeout, retries, and optional circuit breaker
    """
    
    # Add timeout to kwargs if not present
    if 'timeout' not in kwargs:
        kwargs['timeout'] = timeout
    
    last_exception = None
    
    for attempt in range(max_retries + 1):
        try:
            if circuit_breaker_name:
                breaker = get_circuit_breaker(circuit_breaker_name)
                return breaker.call(func, *args, **kwargs)
            else:
                return func(*args, **kwargs)
                
        except Exception as e:
            last_exception = e
            
            if attempt < max_retries:
                wait_time = (attempt + 1) * 2  # Exponential backoff
                logger.warning(f"API call failed (attempt {attempt + 1}/{max_retries + 1}), retrying in {wait_time}s: {str(e)}")
                time.sleep(wait_time)
            else:
                logger.error(f"API call failed after {max_retries + 1} attempts: {str(e)}")
    
    # If we get here, all retries failed
    raise last_exception

# Convenience functions for common timeouts
def embl_api_call(func: Callable, *args, **kwargs) -> Any:
    """Make EMBL API call with appropriate timeout and circuit breaker"""
    return safe_api_call(
        func, *args,
        timeout=15,
        max_retries=2,
        circuit_breaker_name="embl_api",
        **kwargs
    )

def uniprot_api_call(func: Callable, *args, **kwargs) -> Any:
    """Make UniProt API call with appropriate timeout and circuit breaker"""
    return safe_api_call(
        func, *args,
        timeout=20,
        max_retries=2, 
        circuit_breaker_name="uniprot_api",
        **kwargs
    )

def blast_call(func: Callable, *args, **kwargs) -> Any:
    """Make BLAST call with longer timeout"""
    return safe_api_call(
        func, *args,
        timeout=120,  # BLAST can take longer
        max_retries=1,
        circuit_breaker_name="blast_api",
        **kwargs
    )