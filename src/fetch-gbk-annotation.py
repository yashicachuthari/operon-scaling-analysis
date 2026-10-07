#!/usr/bin/env python

'''
fetch-gbk-annotation.py by Rohan Maddamsetti and Irida Shyti.

This script reads in ../results/complete-prokaryotes-with-plasmids.txt.

NOTE: for path names to be processed properly, this script must be run
from the src/ directory as python fetch-gbk-annotation.py.
'''


import urllib.request
from os.path import basename, exists
import os
import subprocess
import sys
from tqdm import tqdm

## imports for parallelization.
import time
import logging
import asyncio
import shutil
from contextlib import asynccontextmanager
##from tqdm.asyncio import tqdm as async_tqdm

class RateLimiter:
    """Rate limiter to prevent overwhelming NCBI servers."""
    
    def __init__(self, calls_per_minute=10):
        self.calls_per_minute = calls_per_minute
        self.interval = 60.0 / calls_per_minute  # seconds between calls
        self.last_call_time = 0
        self.lock = asyncio.Lock()
    
    @asynccontextmanager
    async def limit(self):
        """Context manager to limit the rate of API calls."""
        async with self.lock:
            # Calculate time since last call
            current_time = time.time()
            time_since_last_call = current_time - self.last_call_time
            
            # If we need to wait to respect rate limit
            if time_since_last_call < self.interval:
                wait_time = self.interval - time_since_last_call
                logging.debug(f"Rate limiting: waiting {wait_time:.2f} seconds")
                await asyncio.sleep(wait_time)
            
            # Update last call time
            self.last_call_time = time.time()
        
        try:
            # Allow the caller to proceed
            yield
        finally:
            pass  # No cleanup needed


def configure_logging(log_file):
    try:
        logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(log_file),
            logging.StreamHandler()  # Keep console output too
        ]
    )
    except IOError:
        ## If we can't write to the log file, just log to console
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.StreamHandler()
            ]
        )
        print("WARNING: Could not write to log file. Logging to console only.")

    logging.info("Starting the pipeline...")
    return

        
def create_refseq_accession_to_ftp_path_tuples(prokaryotes_with_plasmids_file):
    refseq_accession_to_ftp_path_tuples = list()
    with open(prokaryotes_with_plasmids_file, "r") as prok_with_plasmids_file_obj:
        for i, line in enumerate(prok_with_plasmids_file_obj):
            if i == 0: continue  # skip the header
            
            fields = line.strip().split("\t")                
            ## Get the accession field (5th from end) and turn GCA Genbank IDs into GCF RefSeq IDs
            refseq_id = fields[-5].replace("GCA", "GCF")
            
            ## Get the ftp_url field (3rd from end) and make sure we turn the GCA Genbank URL
            ## into the GCF RefSeq FTP URL
            ftp_url = fields[-3].replace("GCA", "GCF")

            ## Check for valid IDs and URLs (suppressed RefSeq IDs have a '-' as a blank placeholder)
            if refseq_id.startswith("GCF") and refseq_id in ftp_url:
                my_tuple = (refseq_id, ftp_url)
                refseq_accession_to_ftp_path_tuples.append(my_tuple)
                
    return refseq_accession_to_ftp_path_tuples


def reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
    with open(md5_file, "r") as checksum_fh:
        target_string = "_genomic.gbff.gz"
        for line in checksum_fh:
            if target_string in line:          
                my_target_checksum, my_target_filename = line.strip().split()
                break
    if sys.platform == "darwin":
        my_md5_cmd = "md5" ## run md5 on my mac,
    elif sys.platform == "linux":
        my_md5_cmd = "md5sum" ## but run md5sum on DCC (linux)
    else:
        raise AssertionError("UNKNOWN PLATFORM")

    ## run md5 on the local file and get the output.
    md5_call = subprocess.run([my_md5_cmd, gbff_gz_file], capture_output=True, text=True)

    if sys.platform == "darwin": ## then the checksum is the LAST 'word' in the output.
        my_md5_checksum = md5_call.stdout.split()[-1].strip()
    elif sys.platform == "linux": ## then the checksum is the FIRST 'word' in the output.
        my_md5_checksum = md5_call.stdout.split()[0].strip()
    else:
        raise AssertionError("UNKNOWN PLATFORM")

    ## verify that the checksums match.
    return my_md5_checksum == my_target_checksum


async def download_single_genome(ftp_path, reference_genome_dir):
    """Download a single genome and its MD5 file"""
    my_full_accession = basename(ftp_path)
    my_base_filename = my_full_accession + "_genomic.gbff.gz"
    
    # Files on the NCBI FTP site to download
    gbff_ftp_path = os.path.join(ftp_path, my_base_filename)
    md5_ftp_path = os.path.join(ftp_path, "md5checksums.txt")
    
    # Local paths
    gbff_gz_file = os.path.join(reference_genome_dir, my_base_filename)
    md5_file = os.path.join(reference_genome_dir, my_full_accession + "_md5checksums.txt")

    # Check if files exist and are valid
    if exists(gbff_gz_file) and exists(md5_file):
        if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
            print(f"{gbff_gz_file} SUCCEEDED.")
            return True
        else:
            os.remove(gbff_gz_file)
            os.remove(md5_file)

    # Try downloading up to 5 times
    attempts = 5
    while attempts > 0:
        try:
            await asyncio.wait_for(asyncio.to_thread(urllib.request.urlretrieve, gbff_ftp_path, filename=gbff_gz_file), timeout=120)
            await asyncio.wait_for(asyncio.to_thread(urllib.request.urlretrieve, md5_ftp_path, filename=md5_file), timeout=120)
            
            if reference_genome_passes_md5_checksum(gbff_gz_file, md5_file):
                print(f"{gbff_gz_file} SUCCEEDED.")
                return True
            else:
                if exists(gbff_gz_file):
                    os.remove(gbff_gz_file)
                if exists(md5_file):
                    os.remove(md5_file)
                
        except Exception as e:
            print(f"Attempt {6-attempts} failed for {gbff_gz_file}: {str(e)}")
            if exists(gbff_gz_file):
                os.remove(gbff_gz_file)
            if exists(md5_file):
                os.remove(md5_file)
        
            if "550" in str(e) or "No such file or directory" in str(e):
                break

        attempts -= 1
    
    print(f"{gbff_gz_file} FAILED after all attempts")
    return False


async def async_download(ftp_paths, reference_genome_dir, max_concurrent=10):
    """Download genomes in parallel using asyncio task pool with rate limiting"""
    ## Create a semaphore to limit concurrent downloads
    ## Adjust based on your network capacity
    semaphore = asyncio.Semaphore(max_concurrent)
    
    # Create a rate limiter for NCBI
    rate_limiter = RateLimiter(calls_per_minute=20)  # Adjust as needed
    
    async def download_with_limits(ftp_path):
        """Download a single genome with rate and concurrency limits"""
        async with semaphore:  # Limit concurrent downloads
            async with rate_limiter.limit():  # Respect rate limits
                return await download_single_genome(ftp_path, reference_genome_dir)
    
    # Create tasks for all downloads
    tasks = []
    for ftp_path in ftp_paths:
        task = asyncio.create_task(download_with_limits(ftp_path))
        tasks.append(task)
    
    # Use tqdm to show progress
    results = []
    for task in tqdm(asyncio.as_completed(tasks), total=len(tasks)):
        result = await task
        results.append(result)
    
    return results


def fetch_gbk_annotation(refseq_accession_to_ftp_path_tuples, reference_genome_dir):
    """Download reference genomes for each genome in the RunID table"""
    ## Create reference genome directory if it doesn't exist
    os.makedirs(reference_genome_dir, exist_ok=True)
    
    ## Look up the FTP URLs for each refseq id
    ## IMPORTANT: we have to check to see if the ftp path exists;
    ## this is not true for suppressed entries in RefSeq,
    ## which have an ftp_path == '-'
    ftp_paths = []
    for refseq_id, ftp_path in refseq_accession_to_ftp_path_tuples:
        acc = basename(ftp_path)
        gbff = os.path.join(reference_genome_dir, acc + "_genomic.gbff.gz")
        md5f = os.path.join(reference_genome_dir, acc + "_md5checksums.txt")

        if exists(gbff):
            continue

        ftp_paths.append(ftp_path)

    print(f"Remaining genomes to process: {len(ftp_paths)}")
    ## Run the async download
    asyncio.run(async_download(ftp_paths, reference_genome_dir))
    return


def main():

    ## Configure logging
    log_file = f"../results/gbk_annotation_download.log"
    configure_logging(log_file)

    prokaryotes_with_plasmids_file = "../results/complete-prokaryotes-with-plasmids.txt"
    refseq_accession_to_ftp_path_tuples = create_refseq_accession_to_ftp_path_tuples(prokaryotes_with_plasmids_file)

    reference_genome_dir = "../results/gbk-annotation/"
    fetch_gbk_annotation(refseq_accession_to_ftp_path_tuples, reference_genome_dir)
    return


if __name__ == "__main__":
    main()

