#!/usr/bin/env bash
root_dir=/ROOT_DIR/sanjar/Projects/overbinders/data/DL_analysis/dl_runs/

sbatch --time=100:00:00 \
   --partition=gpu \
   --gres=gpu:p100:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/HepG2_1000_14/std.err \
   --output=$root_dir/HepG2_1000_14/std.out \
   --job-name=H_1k_14 \
   ROOT_DIR/keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py HepG2 1000 14


sbatch --time=100:00:00 \
   --partition=gpu \
   --gres=gpu:p100:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/K562_1000_14/std.err \
   --output=$root_dir/K562_1000_14/std.out \
   --job-name=K_1K_14 \
   ROOT_DIR/keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py K562 1000 14


sbatch --time=50:00:00 \
   --partition=gpu \
   --gres=gpu:p100:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/HepG2_400_14/std.err \
   --output=$root_dir/HepG2_400_14/std.out \
   --job-name=H_400_14 \
   ROOT_DIR/keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py HepG2 400 14

sbatch --time=50:00:00 \
   --partition=gpu \
   --gres=gpu:p100:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/K562_400_14/std.err \
   --output=$root_dir/K562_400_14/std.out \
   --job-name=K_400_14 \
   ROOT_DIR/keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py K562 400 14

exit


sbatch --time=20:00:00 \
   --partition=gpu \
   --gres=gpu:v100x:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/HepG2_400/std.err \
   --output=$root_dir/HepG2_400/std.out \
   --job-name=HepG2_400 \
   keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py HepG2 400 3

sbatch --time=20:00:00 \
   --partition=gpu \
   --gres=gpu:v100x:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/K562_400/std.err \
   --output=$root_dir/K562_400/std.out \
   --job-name=K562_400 \
   keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py K562 400 3

sbatch --time=20:00:00 \
   --partition=gpu \
   --gres=gpu:v100x:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/HepG2_1000/std.err \
   --output=$root_dir/HepG2_1000/std.out \
   --job-name=HepG2_1k \
   keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py HepG2 1000 3

sbatch --time=20:00:00 \
   --partition=gpu \
   --gres=gpu:v100x:1 \
   --cpus-per-task=16 \
   --mem=50g \
   --error=$root_dir/K562_1000/std.err \
   --output=$root_dir/K562_1000/std.out \
   --job-name=K562_1k \
   keras python ROOT_DIR/Projects/enhancers/code/overbinders/dl_analysis/train_model.py K562 1000 3





