#!/bin/bash
#
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=5G
#SBATCH --partition=slow
#SBATCH --time=05:00

# 10x10
srun ./lcs_distributed -x selvyjg0he -y lx3oeh8yso

# 25x25
srun ./lcs_distributed -x qqz96admjoiekgvywb7n55qr0 -y u8dn8wkxiovwa382m77xp137m

# 50x50
srun ./lcs_distributed \
    -x 82t43c8ndm1lbsebl4jfuutw224s9ubzi78k65tidd7b2p9l10 \
    -y ikdugd5uq3rox5lwz63adyvihu2cuo4dkvi1m7p9f7w1nk4za0 
