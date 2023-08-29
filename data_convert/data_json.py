"""
从 data/2011_data_one和 2011_genia_parse_one_split_bio中读取数据
经处理后，写入到 data/2011_json_data中
"""
import os
from process_data import delete_nostandard as dns
import data_utils as data_utils
import json_utils as json_utils
import json

'''
两种数据处理方式: 两个输出序列，第一个只包含实体的输出序列， 第二个包含论元的输出

读原始数据，读genia数据
(1)删除掉3个类型的trigger(2)删除重复事件(4)删除跨句子事件
(3)保留多词组成trigger引发的事件, 这是与gcn_input.py的不同之处
data_path 是一个句子的语料,用于识别trigger和join protein
data_two_path的目的是为了识别trigger与argument的关系以及binding argument
注意 protein 使用的标签系统是BIO, trigger使用的标签系统也是BIO(2020年11月09日)
'''
def read_data(input_path, genia_parse_path, output_path):
    files = os.listdir(input_path)
    all_event_num = 0
    invalid_event_num = 0
    valid_event_num = 0
    sent_json_form = {}
    with open(output_path, "w") as write_fp:
        for file_name in files:
            with open(input_path + "/" + file_name, "r") as data_fp, open(genia_parse_path + "/" + file_name, "r") as genia_parse_fp:
                old_sent = ""
                sent_idx = 0
                str_prot_set = []
                str_prot_idxs = []
                str_jprot_set = []
                str_trig_set = []
                str_event_set = []
                str_event_idx = []

                for line in data_fp:
                    line = line.strip()

                    # if file_name == "PMID-7505113.txt" and line.startswith("Tax co-transfected with reporter constructs"):
                    #     print("#############################")
                    if not len(line):
                        # 处理genia_data中的数据
                        gword_set = []
                        gpos_set = []
                        gprot_set = []
                        gjprot_set = []
                        cjprot_set = []
                        ctrig_set = []
                        trigtype_set = []

                        genia_parse_line = genia_parse_fp.readline()
                        genia_parse_line = genia_parse_line.strip().strip('\n')
                        while len(genia_parse_line) != 0:
                            genia_parse = genia_parse_line.split("\t###\t")
                            # genia_line信息处理，为识别trigger做准备
                            # (word+"\t"+pos+"\t"+fea1+"\t"+fea2+"\t"+prot+"\t"+cjprot+"\t"+ctrig+"\t"+jprot+"\t"+trigtype+"\n")
                            genia_line = genia_parse[0]
                            genia_info = genia_line.split("\t")
                            gword_set.append(genia_info[0])
                            gpos_set.append(genia_info[1])
                            gprot_set.append(genia_info[4])
                            cjprot_set.append(genia_info[5])
                            ctrig_set.append(genia_info[6])
                            gjprot_set.append(genia_info[7])
                            if genia_info[8] == "Entity":
                                trigtype_set.append('Other')
                            else:
                                trigtype_set.append(genia_info[8])

                            genia_parse_line = genia_parse_fp.readline()
                            genia_parse_line = genia_parse_line.strip().strip('\n')

                        new_sent = " ".join(gword_set)
                        assert old_sent.replace(" ", "") == new_sent.replace(" ", "")

                        # --------(1)删除掉四个类型的trigger(2)删除重复事件(3)删除跨句子事件(3)删除多词trigger引发的事件--------------
                        valid_event_set = dns.delete_invalid_event_two(str_event_set, str_prot_idxs, str_event_idx)
                        # -----（4）一个trigger词有多个类型， 保留一个类型(保留简单类型)
                        valid_trig_set, valid_event_set = dns.delete_invalid_event(str_trig_set, valid_event_set)
                        # ----------------------------------------------------------
                        invalid_event_num += len(str_event_set) - len(valid_event_set)
                        valid_event_num += len(valid_event_set)
                        # 将str_prot, str_trig, str_event转换成event
                        prot_dict = data_utils.prot_str2entity(old_sent, new_sent, str_prot_set)
                        #jprot_dict = utils.prot_str2entity(old_sent, new_sent, str_jprot_set)
                        trig_dict = data_utils.trig_str2entity(old_sent, new_sent, valid_trig_set)
                        event_dict = data_utils.event_str2entity(valid_event_set, trig_dict, prot_dict)
                        assert len(valid_event_set) == len(event_dict)

                        sent_json_form["doc_id"] = file_name
                        sent_json_form["sent_id"] = "S"+str(sent_idx)
                        sent_json_form["tokens"] = gword_set
                        sent_json_form["sentence"] = new_sent

                        # 添加 entity_mentions and event_mentions
                        json_utils.data2json(prot_dict, trig_dict, event_dict, sent_json_form)

                        write_fp.write(json.dumps(sent_json_form) + "\n")

                        old_sent = ""
                        str_prot_set.clear()
                        str_jprot_set.clear()
                        str_trig_set.clear()
                        str_event_set.clear()
                        str_event_idx.clear()
                        str_prot_idxs.clear()

                    elif line.startswith("#"):
                        str_prot_set.append(line)
                        prot_info = line.split()
                        str_prot_idxs.append(prot_info[1])

                    elif line.startswith("$"):
                        str_jprot_set.append(line)

                    elif line.startswith("@"):# @ T20 Regulation S2 160 172 306 318 deregulation
                        trig_info = line.split()
                        if trig_info[2] not in ["Entity", 'Anaphora', 'Ubiquitination', 'Protein_modification']:
                            str_trig_set.append(line)

                    elif line.startswith("%"):
                        all_event_num += 1
                        str_event, event_idx = data_utils.process_strevent(line)
                        if str_event not in str_event_set: ## 删除重复event
                            str_event_set.append(str_event)
                            str_event_idx.append(event_idx)
                    else:
                        old_sent = line
                        sent_idx += 1

    #comp_prev_next_sent(file_set)
    print("所有事件个数:", all_event_num)
    print("跨句子事件个数：", invalid_event_num)
    print("有效事件个数：", valid_event_num)
    print("#################################")

'''
读取test集中数据
'''
def read_test_data(data_path, parse_path, output_path):
    files = os.listdir(data_path)
    sent_json_form = {}
    with open(output_path, "w") as output_fp:
        for file_name in files:
            with open(data_path + "/" + file_name, "r") as data_fp, open(parse_path + "/" + file_name, "r") as parse_fp:
                old_sent = ""
                sent_idx = 0
                str_prot_set = list()

                for line in data_fp:
                    line = line.strip()
                    if not len(line):
                        ## genia 信息数据
                        gword_set = list()
                        gpos_set = list()
                        gprot_set = list()

                        genia_parse_line = parse_fp.readline()
                        genia_parse_line = genia_parse_line.strip().strip('\n')
                        while len(genia_parse_line) != 0:
                            genia_parse = genia_parse_line.split("\t###\t")
                            # genia_line: Down-regulation	NN	B-NP	O	Other
                            # (word+"\t"+pos+"\t"+fea1+"\t"+fea2+"\t"+prot")
                            genia_line = genia_parse[0].strip()
                            genia_info = genia_line.split("\t")
                            gword_set.append(genia_info[0])
                            gpos_set.append(genia_info[1])
                            gprot_set.append(genia_info[4])

                            genia_parse_line = parse_fp.readline()
                            genia_parse_line = genia_parse_line.strip().strip('\n')

                        new_sent = " ".join(gword_set)
                        if old_sent.replace(" ", "") != new_sent.replace(" ", ""):
                            print("#########################################")
                            print(old_sent.replace(" ", ""))
                            print(new_sent.replace(" ", ""))
                            print("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
                        assert old_sent.replace(" ", "") == new_sent.replace(" ", "")

                        prot_dict = data_utils.prot_str2entity(old_sent, new_sent, str_prot_set)

                        sent_json_form["doc_id"] = file_name
                        sent_json_form["sent_id"] = "S"+str(sent_idx)
                        sent_json_form["tokens"] = gword_set
                        sent_json_form["sentence"] = new_sent
                        # 添加 entity_mentions and event_mentions
                        json_utils.testdata2json(prot_dict, sent_json_form)
                        output_fp.write(json.dumps(sent_json_form) + "\n")

                        old_sent = ""
                        str_prot_set.clear()

                    elif line.startswith("#"):
                        str_prot_set.append(line)

                    elif line.startswith("$") or line.startswith("@") or line.startswith("%"):
                        continue
                    else:
                        old_sent = line
                        sent_idx += 1


if __name__ == '__main__':
    absolute_path = "//"

    train_data = absolute_path + "data/2011_data_one/2011_train"
    train_genia = absolute_path + "data/2011_genia_parse_one_split_bio/2011_train"
    train_processed = absolute_path + "data/2011_json_data/train.json" # 转换为增广后的数据
    read_data(train_data, train_genia, train_processed)

    devel_data = absolute_path + "data/2011_data_one/2011_devel"
    devel_genia = absolute_path + "data/2011_genia_parse_one_split_bio/2011_devel"
    devel_processed = absolute_path + "data/2011_json_data/devel.json"
    read_data(devel_data, devel_genia, devel_processed)

    test_data = absolute_path + "data/2011_data_one/2011_test"
    test_genia = absolute_path + "data/2011_genia_parse_one_split_bio/2011_test"
    test_processed = absolute_path + "data/2011_json_data/test.json"
    read_test_data(test_data, test_genia, test_processed)
