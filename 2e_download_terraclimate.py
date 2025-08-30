#!/usr/bin/env python3
"""
Download TerraClimate data for 2001-2019
Variables: precipitation, soil moisture, temperature, VPD
Downloads to: input_data/terraclimate_raw/
"""

import os
import requests
import yaml
from time import sleep

def load_config():
    """Load configuration from config.yaml"""
    try:
        with open('config.yaml', 'r') as file:
            return yaml.safe_load(file)
    except:
        return {}

def download_file(url, filename):
    """Download with progress bar and better error handling"""
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36'
    }
    
    try:
        response = requests.get(url, stream=True, headers=headers, timeout=30)
        response.raise_for_status()
        
        total_size = int(response.headers.get('content-length', 0))
        
        # Check if we're getting an HTML error page
        content_type = response.headers.get('content-type', '')
        if 'text/html' in content_type:
            print(" ✗ Error: Received HTML instead of data")
            return False
        
        with open(filename, 'wb') as file:
            downloaded = 0
            for chunk in response.iter_content(chunk_size=8192):
                file.write(chunk)
                downloaded += len(chunk)
                if total_size > 0:
                    percent = (downloaded / total_size) * 100
                    print(f"\r  {percent:.1f}% complete", end='')
        print()
        return True
        
    except requests.exceptions.RequestException as e:
        print(f" ✗ Error: {e}")
        return False

def check_file_is_netcdf(filepath):
    """Quick check if file is actually NetCDF"""
    try:
        with open(filepath, 'rb') as f:
            header = f.read(4)
            # NetCDF files start with 'CDF\x01' or 'CDF\x02' or '\x89HDF'
            return header.startswith(b'CDF') or header.startswith(b'\x89HDF')
    except:
        return False

def main():
    # Define variables to download
    variables = ['ppt', 'soil', 'tmax', 'vpd']
    start_year = 2001
    end_year = 2019
    
    # Create download directory
    download_dir = "input_data/terraclimate_raw"
    os.makedirs(download_dir, exist_ok=True)
    
    # Clean up any 4KB files from previous attempt
    print("Cleaning up invalid files from previous attempt...")
    cleaned = 0
    for f in os.listdir(download_dir):
        if f.endswith('.nc'):
            filepath = os.path.join(download_dir, f)
            if os.path.getsize(filepath) < 1000000:  # Less than 1MB
                os.remove(filepath)
                cleaned += 1
    if cleaned > 0:
        print(f"Removed {cleaned} invalid files\n")
    
    # The correct URL pattern (found from TerraClimate website)
    base_url = "https://climate.northwestknowledge.net/TERRACLIMATE-DATA"
    
    print(f"Starting TerraClimate download ({start_year}-{end_year})")
    print(f"Variables: {', '.join(variables)}")
    print(f"Download directory: {download_dir}")
    print(f"Total files to download: {len(variables) * (end_year - start_year + 1)}")
    print(f"Using URL pattern: {base_url}/TerraClimate_VARIABLE_YEAR.nc")
    print()
    
    # Track overall progress
    total_files = len(variables) * (end_year - start_year + 1)
    completed_files = 0
    failed_files = []
    
    for var in variables:
        print(f"\nDownloading {var} data...")
        for year in range(start_year, end_year + 1):
            completed_files += 1
            filename = f"TerraClimate_{var}_{year}.nc"
            filepath = os.path.join(download_dir, filename)
            
            # Check if valid file already exists
            if os.path.exists(filepath):
                size = os.path.getsize(filepath) / 1e6
                if size > 10 and check_file_is_netcdf(filepath):  # Expect files > 10MB
                    print(f"  [{completed_files}/{total_files}] {year}: already exists ({size:.1f} MB)")
                    continue
                else:
                    print(f"  [{completed_files}/{total_files}] {year}: invalid file, re-downloading...")
                    os.remove(filepath)
            
            url = f"{base_url}/{filename}"
            print(f"  [{completed_files}/{total_files}] {year}: downloading...", end='')
            
            success = download_file(url, filepath)
            
            if success and os.path.exists(filepath):
                size = os.path.getsize(filepath) / 1e6
                if size > 10 and check_file_is_netcdf(filepath):
                    print(f" ✓ ({size:.1f} MB)")
                else:
                    print(f" ✗ Invalid file ({size:.1f} MB)")
                    os.remove(filepath)
                    failed_files.append(filename)
            else:
                if os.path.exists(filepath):
                    os.remove(filepath)
                failed_files.append(filename)
            
            # Small delay to be nice to the server
            sleep(0.5)
    
    # Summary
    print("\nDownload summary:")
    for var in variables:
        files_found = len([f for f in os.listdir(download_dir) 
                          if f.startswith(f"TerraClimate_{var}_") and f.endswith('.nc')])
        print(f"  {var}: {files_found}/{end_year - start_year + 1} files")
    
    if any(f.endswith('.nc') for f in os.listdir(download_dir)):
        total_size = sum(os.path.getsize(os.path.join(download_dir, f)) 
                         for f in os.listdir(download_dir) if f.endswith('.nc'))
        print(f"\nTotal size: {total_size/1e9:.2f} GB")
    
    if failed_files:
        print(f"\nFailed downloads ({len(failed_files)} files):")
        for f in failed_files[:10]:  # Show first 10
            print(f"  - {f}")
        if len(failed_files) > 10:
            print(f"  ... and {len(failed_files) - 10} more")
    
    print("\nDownload completed!")
    if not failed_files:
        print("Next step: Run 2f_process_terraclimate.R to create 2001-2019 averages")

if __name__ == "__main__":
    main()