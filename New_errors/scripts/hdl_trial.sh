#PBS -lwalltime=10:00:00
#PBS -lselect=1:ncpus=4:mem=32gb

module load anaconda3/personal

cd ~/Summer_project/Scripts

Rscript hdl_trial_server.R
