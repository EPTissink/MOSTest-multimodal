from string import Template
import subprocess
import os

ref_file = "/path/ref9545380_1kgPhase3eur_LDr2p1.mat"
trait_folder = "/path/sumstats"

trait_file1 = ['PGC_MDD_2018_no23andMe_noUKBB.mat', 'PGC_ADHD_2017_EUR.mat', 'PGC_ASD_2017_iPSYCH.mat']
trait_name1 = ['MDD', 'ADHD', 'ASD']
trait_file2 = ['multimodal_eig26.most.orig.mat']
trait_name2 = ['MOSTest-multimodal']

output_dir_tmp = "/cluster/projects/p33/users/alexeas/elleke/pleioprs/condfdr_results"

stat_type = "condfdr"
fdr_thresh = "0.05"
randprune_n = "100"
exclude_chr_pos = "6 25000000 33000000; 8 7200000 12500000"

config_template_file = "config_template.txt"
with open(config_template_file) as f:
    config_template = Template(f.read())

par_dict = dict(ref_file=ref_file, trait_folder=trait_folder, stat_type=stat_type,
                fdr_thresh=fdr_thresh, randprune_n=randprune_n, exclude_chr_pos=exclude_chr_pos)

for tf1, tn1  in zip(trait_file1, trait_name1):
    for tf2, tn2 in zip(trait_file2, trait_name2):
        print(f'Running {stat_type} for {tn1} and {tn2}')
        output_dir = os.path.join(output_dir_tmp, f'{stat_type}_{tn1}_{tn2}')
        par_dict["out_dir"] = output_dir
        par_dict["trait_file1"] = tf1
        par_dict["trait_name1"] = tn1
        par_dict["trait_file2"] = tf2
        par_dict["trait_name2"] = tn2
        config = config_template.substitute(par_dict)
        with open('config.txt', 'w') as f:
            f.write(config)
        command = ['matlab', '-nodisplay', '-nosplash', '<', 'runme.m']
        subprocess.run(command)

print('Done')

