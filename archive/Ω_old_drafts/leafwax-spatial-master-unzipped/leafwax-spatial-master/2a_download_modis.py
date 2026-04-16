#!/usr/bin/env python3
"""
MODIS MCD12Q1 Data Download Script
Downloads MCD12Q1 Collection 6.1 data from 2001-2019 using NASA earthaccess
"""

import earthaccess
import os
import yaml

def load_config():
    """Load configuration from config.yaml"""
    with open('config.yaml', 'r') as file:
        return yaml.safe_load(file)
        
def main():
    # Load configuration
    config = load_config()
    modis_config = config['modis']
    
    print(f"Starting MODIS {modis_config['product']} download ({modis_config['start_year']}-{modis_config['end_year']})")
    
    # Create download directory from config
    os.makedirs(modis_config['download_dir'], exist_ok=True)
    
    # Authenticate with NASA Earthdata
    print("Authenticating with NASA Earthdata...")
    auth = earthaccess.login()
    
    if not auth:
        print("Authentication failed! Please check your NASA Earthdata credentials.")
        return
    
    print("✓ Authentication successful")
    
    # Search for MCD12Q1 data using config parameters
    print(f"Searching for {modis_config['product']} data ({modis_config['start_year']}-{modis_config['end_year']})...")
    results = earthaccess.search_data(
        short_name=modis_config['product'],
        version=modis_config['version'],
        temporal=(f"{modis_config['start_year']}-01-01", f"{modis_config['end_year']}-12-31")
    )
    
    print(f"Found {len(results)} files")
    
    # Download files
    print("Starting download...")
    earthaccess.download(results, modis_config['download_dir'])
    
    print("Download completed!")

if __name__ == "__main__":
    main()