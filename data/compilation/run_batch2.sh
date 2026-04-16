#!/bin/bash
# Sequential orchestrator for confounding v2 batch 2
# Runs rho05 then empirical, shuts down when done

cd /home/shared/leafwax_spatial/spatial_leafwax_model
LOG=logs/batch_v2_orchestrator.log

echo "Sequential batch 2 (treedepth=12) started at $(date)" >> "$LOG"

echo "Starting rho05..." >> "$LOG"
Rscript run_confounding_test_v2.R rho05 > logs/confounding_v2_rho05.log 2>&1
echo "rho05 done at $(date)" >> "$LOG"

echo "Starting empirical..." >> "$LOG"
Rscript run_confounding_test_v2.R empirical > logs/confounding_v2_empirical.log 2>&1
echo "empirical done at $(date)" >> "$LOG"

echo "All v2 batch 2 complete at $(date). Shutting down." >> "$LOG"
sudo shutdown -h now
