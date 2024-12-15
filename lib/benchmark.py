import gc
import time
import psutil
import functools
from typing import Callable, Any
import tracemalloc

def benchmark(func: Callable) -> Callable:
    """
    A decorator that measures execution time and memory usage of a function.
    
    Args:
        func: The function to benchmark
        
    Returns:
        Wrapped function that prints benchmark results after execution
    """
    @functools.wraps(func)
    def wrapper(*args, **kwargs) -> Any:
        # Disable automatic garbage collection
        gc.disable()
        
        # Force garbage collection before starting
        gc.collect()
        
        # Start tracemalloc for Python memory tracking
        tracemalloc.start()
        
        # Get initial process memory (RSS)
        process = psutil.Process()
        process.memory_info()  # First call to ensure cache is warm
        rss_before = process.memory_info().rss
        
        # Start timing
        start_time = time.perf_counter()
        
        # Execute function
        result = func(*args, **kwargs)
        
        # End timing
        end_time = time.perf_counter()
        execution_time = end_time - start_time
        
        # Force garbage collection before final measurement
        gc.collect()
        
        # Get memory measurements
        rss_after = process.memory_info().rss
        current, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        
        # Re-enable automatic garbage collection
        gc.enable()
        
        # Convert all measurements to MB for consistency
        rss_before_mb = rss_before / 1024 / 1024  # MB
        rss_after_mb = rss_after / 1024 / 1024    # MB
        rss_diff = (rss_after_mb - rss_before_mb) # MB
        current_mb = current / 1024 / 1024  # MB
        peak_mb = peak / 1024 / 1024  # MB

        return {
            "ResultObject": result,
            "Time": f"{execution_time:.4f}",
            "RSS": f"{rss_diff:.2f}",
            "PeakMemory": f"{peak_mb:.2f}",
            "FinalMemory": f"{current_mb:.2f}"
        }
    return wrapper