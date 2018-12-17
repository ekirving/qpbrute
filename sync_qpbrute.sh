#!/usr/bin/env bash

rsync -avz --partial \
           --progress \
           --stats \
           --exclude '.git*' \
           --exclude '*.sh' \
           --exclude '*.txt' \
           --exclude '*.py*' \
           --exclude '*.r' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '*tmp*' \
           --exclude '*.gz' \
           --exclude '*.bak' \
           --exclude 'nohup*' \
           --exclude 'graphs/*.log' \
           --exclude 'bayes/*-chain.csv' \
           evan@palaeoprime:~/qpbrute/ ~/Dropbox/Code/qpbrute/
