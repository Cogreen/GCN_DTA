# Date: 2020.11.27

import os
import random
import numpy as np
import json
from collections import OrderedDict
from datetime import datetime
import csv 

def seq_format(proteins_dic, output_dir):
    # datasets = ['kiba', 'davis']
    # for dataset in datasets:
    # seq_path = os.path.join('data', dataset, 'seq_original')
    # fpath = os.path.join('data', dataset, 'proteins.txt')
    # proteins = json.load(fpath, object_pairs_hook=OrderedDict)
    for key, value in proteins_dic.items():
        with open(os.path.join(output_dir, key + '.fasta'), 'w') as f:
            f.writelines('>' + key + '\r\n')
            f.writelines(value + '\r\n')


#count = 0 
def HHblitsMSA(bin_path, db_path, input_dir, output_dir):
    count = 0 
    for fas_file in os.listdir(input_dir):
        HHblitsMSA_start = datetime.now()
        HHblitsMSA_start1 = HHblitsMSA_start.strftime('%Y-%m-%d %H:%M:%S')
        file_name = fas_file.split('.fasta')[0]
        process_file = os.path.join(input_dir, fas_file)
        output_file = os.path.join(output_dir, fas_file.split('.fasta')[0] + '.hhr')  # igore
        output_file_a3m = os.path.join(output_dir, fas_file.split('.fasta')[0] + '.a3m')
        if os.path.exists(output_file) and os.path.exists(output_file_a3m):
            print(output_file, output_file_a3m, 'exist.')
            count += 1
            continue
        print(process_file)
        process_file = process_file.replace('(', '\(').replace(')', '\)')
        print(process_file)
        output_file = output_file.replace('(', '\(').replace(')', '\)')
        print(output_file)
        output_file_a3m = output_file_a3m.replace('(', '\(').replace(')', '\)')
        print(output_file_a3m)
        cmd = bin_path + ' -maxfilt 100000 -realign_max 100000 -d ' + db_path + ' -all -B 100000 -Z 100000 -n 3 -e 0.001 -i ' + process_file + ' -o ' + output_file + ' -oa3m ' + output_file_a3m + ' -cpu 8'
        print(cmd)
        os.system(cmd)
        HHblitsMSA_end = datetime.now()
        HHblitsMSA_end1 = HHblitsMSA_end.strftime('%Y-%m-%d %H:%M:%S')
        HHblitsMSA_duration = HHblitsMSA_end - HHblitsMSA_start
        HHblitsMSA_duration1 = str(HHblitsMSA_duration)
        
        with open("scripts_duration.csv", mode = "a") as f1:
            f1 = csv.writer(f1, delimiter = ',')
            f1.writerow(['HHblitsMSA', file_name, HHblitsMSA_start1, HHblitsMSA_end1, HHblitsMSA_duration1])
        print(count)


def HHfilter(bin_path, input_dir, output_dir):
    file_prefix = []
    # print(input_dir)
    for file in os.listdir(input_dir):
        if 'a3m' not in file:
            continue
        temp_prefix = file.split('.a3m')[0]
        if temp_prefix not in file_prefix:
            file_prefix.append(temp_prefix)
    # random.shuffle(file_prefix)
    # print(len(file_prefix))
    # print(file_prefix)
    for msa_file_prefix in file_prefix:
        HHfilter_start = datetime.now()
        HHfilter_start1 = HHfilter_start.strftime('%Y-%m-%d %H:%M:%S') 
        file_name = msa_file_prefix + '.a3m'
        process_file = os.path.join(input_dir, file_name)
        output_file = os.path.join(output_dir, file_name)
        if os.path.exists(output_file):
            continue
        process_file = process_file.replace('(', '\(').replace(')', '\)')
        output_file = output_file.replace('(', '\(').replace(')', '\)')
        cmd = bin_path + ' -id 90 -i ' + process_file + ' -o ' + output_file
        print(cmd)
        os.system(cmd)
        HHfilter_end = datetime.now()
        HHfilter_end1 = HHfilter_end.strftime('%Y-%m-%d %H:%M:%S')
        HHfilter_duration = HHfilter_end - HHfilter_start
        HHfilter_duration1 = str(HHfilter_duration)
        with open("scripts_duration.csv", mode = "a") as f1:
            f1 = csv.writer(f1, delimiter = ',')
            f1.writerow(['HHfilter',msa_file_prefix, HHfilter_start1, HHfilter_end1, HHfilter_duration1])


def reformat(bin_path, input_dir, output_dir):
    # print('reformat')
    for a3m_file in os.listdir(input_dir):
        file_name = a3m_file.split('.a3m')[0]
        reformat_start = datetime.now()
        reformat_start1 = reformat_start.strftime('%Y-%m-%d %H:%M:%S')
        process_file = os.path.join(input_dir, a3m_file)
        output_file = os.path.join(output_dir, a3m_file.split('.a3m')[0] + '.fas')
        if os.path.exists(output_file):
            continue
        process_file = process_file.replace('(', '\(').replace(')', '\)')
        output_file = output_file.replace('(', '\(').replace(')', '\)')
        cmd = bin_path + ' ' + process_file + ' ' + output_file + ' -r'
        print(cmd)
        os.system(cmd)
        reformat_end = datetime.now()
        reformat_end1 = reformat_end.strftime('%Y-%m-%d %H:%M:%S')
        reformat_duration = reformat_end - reformat_start
        reformat_duration1 = str(reformat_duration)
        with open("scripts_duration.csv", mode = "a") as f1:
            f1 = csv.writer(f1, delimiter = ',')
            f1.writerow(['reformat', file_name, reformat_start1, reformat_end1, reformat_duration1])


def convertAlignment(bin_path, input_dir, output_dir):
    # print('convertAlignment')
    for fas_file in os.listdir(input_dir):
        convertAlignment_start = datetime.now()
        file_name = fas_file.split('.fas')[0]
        convertAlignment_start1 = convertAlignment_start.strftime('%Y-%m-%d %H:%M:%S')
        process_file = input_dir + '/' + fas_file
        output_file = output_dir + '/' + fas_file.split('.fas')[0] + '.aln'
        if os.path.exists(output_file):
            continue
        process_file = process_file.replace('(', '\(').replace(')', '\)')
        output_file = output_file.replace('(', '\(').replace(')', '\)')
        cmd = 'python ' + bin_path + ' ' + process_file + ' fasta ' + output_file
        print(cmd)
        os.system(cmd)
        convertAlignment_end = datetime.now()
        convertAlignment_end1 = convertAlignment_end.strftime('%Y-%m-%d %H:%M:%S')
        convertAlignment_duration = convertAlignment_end - convertAlignment_start
        convertAlignment_duration1 = str(convertAlignment_duration)
        with open("scripts_duration.csv", mode = "a") as f1:
            f1 = csv.writer(f1, delimiter = ',')
            f1.writerow(['convertAlignment', file_name, convertAlignment_start1, convertAlignment_end1, convertAlignment_duration1])


def alnFilePrepare():
    import json
    from collections import OrderedDict
    print('aln file prepare ...')
    datasets = ['davis', 'kiba']
    for dataset in datasets:
        seq_dir = os.path.join('data', dataset, 'seq')
        msa_dir = os.path.join('data', dataset, 'msa')
        filter_dir = os.path.join('data', dataset, 'hhfilter')
        reformat_dir = os.path.join('data', dataset, 'reformat')
        aln_dir = os.path.join('data', dataset, 'aln')
        pconsc4_dir = os.path.join('data', dataset, 'pconsc4')
        protein_path = os.path.join('data', dataset)
        proteins = json.load(open(os.path.join(protein_path, 'proteins.txt')), object_pairs_hook=OrderedDict)

        if not os.path.exists(seq_dir):
            os.makedirs(seq_dir)
        if not os.path.exists(msa_dir):
            os.makedirs(msa_dir)
        if not os.path.exists(filter_dir):
            os.makedirs(filter_dir)
        if not os.path.exists(reformat_dir):
            os.makedirs(reformat_dir)
        if not os.path.exists(aln_dir):
            os.makedirs(aln_dir)

    #    HHblits_bin_path = '..../tool/hhsuite/bin/hhblits'  # HHblits bin path
    #    HHblits_db_path = '..../dataset/uniclust/uniclust30_2018_08/uniclust30_2018_08'  # hhblits dataset for msa
    #    HHfilter_bin_path = '..../tool/hhsuite/bin/hhfilter'  # HHfilter bin path
    #    reformat_bin_path = '..../tool/hhsuite/scripts/reformat.pl'  # reformat bin path
    #    convertAlignment_bin_path = '..../tool/CCMpred/scripts/convert_alignment.py'  # ccmpred convertAlignment bin path
        HHblits_bin_path = '/home/ailon/tool/hhsuite/bin/hhblits'  # HHblits bin path
        HHblits_db_path = '/mnt/ailon/dataset/uniclust30_2018_08/uniclust30_2018_08'  # hhblits dataset for msa
        HHfilter_bin_path = '/home/ailon/tool/hhsuite/bin/hhfilter'  # HHfilter bin path
        reformat_bin_path = '/home/ailon/tool/hhsuite/scripts/reformat.pl'  # reformat bin path
        convertAlignment_bin_path = '/home/ailon/tool/CCMpred/scripts/convert_alignment.py'  # ccmpred convertAlignment bin path

        # check the programs used for the script
        if not os.path.exists(HHblits_bin_path):
            raise Exception('Program HHblits was not found. Please specify the run path.')

        if not os.path.exists(HHfilter_bin_path):
            raise Exception('Program HHfilter was not found. Please specify the run path.')

        if not os.path.exists(reformat_bin_path):
            raise Exception('Program reformat was not found. Please specify the run path.')

        if not os.path.exists(convertAlignment_bin_path):
            raise Exception('Program convertAlignment was not found. Please specify the run path.')

        #seq_format(proteins, seq_dir)      
        HHblitsMSA(HHblits_bin_path, HHblits_db_path, seq_dir, msa_dir)
        HHfilter(HHfilter_bin_path, msa_dir, filter_dir)
        reformat(reformat_bin_path, filter_dir, reformat_dir)
        convertAlignment(convertAlignment_bin_path, reformat_dir, aln_dir)

        print('aln file prepare over.')


def pconsc4Prediction():
    import pconsc4
    import keras.backend.tensorflow_backend as KTF
    import tensorflow as tf
    # config = tf.ConfigProto()
    # config.gpu_options.allow_growth = True
    # session = tf.Session(config=config)
    # KTF.set_session(session)
    datasets = ['davis', 'kiba']
    model = pconsc4.get_pconsc4()
    for dataset in datasets:
        aln_dir = os.path.join('data', dataset, 'hhfilter')
        output_dir = os.path.join('data', dataset, 'pconsc4')
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        file_list = os.listdir(aln_dir)
        random.shuffle(file_list)
        inputs = []
        outputs = []
        for file in file_list:
            pconsc4_start_time = datetime.now()
            pconsc4_start_time1 = pconsc4_start_time.strftime('%Y-%m-%d %H:%M:%S')
            input_file = os.path.join(aln_dir, file)
            output_file = os.path.join(output_dir, file.split('.a3m')[0] + '.npy')
            if os.path.exists(output_file):
                # print(output_file, 'exist.')
                continue
            inputs.append(input_file)
            outputs.append(output_file)
            try:
                print('process', input_file)
                pred = pconsc4.predict(model, input_file)
                np.save(output_file, pred['cmap'])
                print(output_file, 'over.')
                pconsc4_end_time = datetime.now()
                pconsc4_end_time1 = pconsc4_end_time.strftime('%Y-%m-%d %H:%M:%S')
                pconsc4_duration = pconsc4_end_time - pconsc4_start_time
                pconsc4_duration1 = str(pconsc4_duration)
                with open("scripts_duration.csv", mode = "a") as f2:                        
                   f2 = csv.writer(f2, delimiter=',')
                   f2.writerow(['PconsC4', file, pconsc_stsart_time1, pconsc_end_time1, pconsc4_duration1])
            except:
                print(output_file, 'error.')




if __name__ == '__main__':
    alnFilePrepare()
    pconsc4Prediction()
