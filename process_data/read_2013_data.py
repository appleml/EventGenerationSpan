import os
import json
import delete_nostandard as dns
import data_utils as util

'''
读原始数据，读parser解析出的数据
event 删除的类型:
(1)删除跨句子事件
(2)删除重复事件 (在读取的时候已经删除啦)
(3) overlap trigger保留一个类型
(4) 删除特殊trigger类型 ['Anaphora', 'Ubiquitination', 'Protein_modification']

Event 保留的类型：
(1) 多词 trigger 引发的事件, (trigger 由多个词组成)
(2) 多类型 trigger 引发的事件， (词由多个trigger类型)

目的：
转换成 
'''
def read_data(data_path, genia_parse_path, output_path, span_output_path):
    files = os.listdir(data_path)
    has_event_sent_num = 0
    changed_event_tree = 0
    sent_json_form = {}
    span_sent_json_form = {}

    for file_name in files:
        filename, extension = os.path.splitext(file_name)
        new_filename = filename + ".json"
        with open(data_path + file_name, "r") as data_fp, open(genia_parse_path + file_name, "r") as genia_parse_fp, \
                open(output_path + new_filename, "w") as writer, open(span_output_path + new_filename, "w") as span_writer:
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
                if not len(line):
                    # 处理genia_data中的数据
                    gword_set = []
                    gpos_set = []
                    gprot_set = []

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

                        genia_parse_line = genia_parse_fp.readline()
                        genia_parse_line = genia_parse_line.strip().strip('\n')

                    new_sent = " ".join(gword_set)
                    if old_sent.replace(" ", "") != new_sent.replace(" ", ""):
                        print(old_sent)
                        print(new_sent)
                        print("$$$$$$$$$$$$$$$$$$$$$$$$$$$44")
                    assert old_sent.replace(" ", "") == new_sent.replace(" ", "")

                    # -----------------------------(1)删除跨句子事件--------------------------------------------
                    valid_event_set = dns.delete_invalid_event_two(str_event_set, str_prot_idxs, str_event_idx)

                    # 将str_prot, str_trig, str_event转换成event
                    prot_dict = util.prot_str2entity(old_sent, new_sent, str_prot_set)
                    # jprot_dict = util.prot_str2entity(old_sent, new_sent, str_jprot_set)
                    trig_dict = util.trig_str2entity(old_sent, new_sent, str_trig_set)
                    event_dict = util.event_str2entity(valid_event_set, trig_dict, prot_dict)
                    assert len(valid_event_set) == len(event_dict)

                    # 生成输入， 带蛋白质的输入
                    augmented_input = util.gen_input_seq(new_sent, prot_dict)
                    # 生成输出， 事件序列
                    event_tree = util.gen_output_seq(file_name, old_sent, trig_dict, event_dict)
                    span_event_tree = util.gen_output_seq_span(file_name, old_sent, trig_dict, event_dict)
                    if len(event_dict) > 0:
                        has_event_sent_num += 1
                    if event_tree != "<extra_id_0> <extra_id_1>":
                        changed_event_tree += 1

                    sent_json_form["text"] = augmented_input
                    sent_json_form["event"] = event_tree
                    writer.write(json.dumps(sent_json_form) + "\n")
                    sent_json_form.clear()

                    # 生成 span
                    span_sent_json_form["text"] = augmented_input
                    span_sent_json_form["event"] = span_event_tree
                    span_writer.write(json.dumps(span_sent_json_form) + "\n")
                    span_sent_json_form.clear()

                    old_sent = ""
                    str_prot_set.clear()
                    str_jprot_set.clear()
                    str_trig_set.clear()
                    str_event_set.clear()
                    str_event_idx.clear()
                    str_prot_idxs.clear()

                elif line.startswith("#") and not (line.startswith("###") or line.startswith("#p<0.05")):
                    str_prot_set.append(line)
                    prot_info = line.split()
                    str_prot_idxs.append(prot_info[1])
                    #assert sent_idx == int(prot_info[3][1:])

                elif line.startswith("$"):
                    str_jprot_set.append(line)

                elif line.startswith("@"):# @ T20 Regulation S2 160 172 306 318 deregulation
                    trig_info = line.split()
                    if trig_info[2] not in ["Entity", 'Anaphora']:
                        str_trig_set.append(line)

                elif line.startswith("%"):
                    str_event, event_idx = util.process_strevent(line)
                    if str_event not in str_event_set: ## 删除重复event
                        str_event_set.append(str_event)
                        str_event_idx.append(event_idx)
                else:
                    old_sent = line
                    sent_idx += 1
                    # print(file_name)
                    # print(old_sent)
                    # print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
                    # if old_sent.startswith("Each bar represents mean+/-SEM of at least three independent"):
                    #     print("@@@@@@@@@@@@@@@@@@@@@@@@")
    #assert changed_event_tree == has_event_sent_num, has_event_sent_num
    print(changed_event_tree, has_event_sent_num)

'''
读取test集中数据
'''
def read_test_data(data_path, parse_path, output_path):
    files = os.listdir(data_path)
    sent_json_form = dict()

    for file_name in files:
        filename, extension = os.path.splitext(file_name)
        new_filename = filename + ".json"
        with open(data_path + file_name, "r") as data_fp, open(parse_path + file_name, "r") as parse_fp, open(output_path+new_filename, "w") as writer:
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
                    assert old_sent.replace(" ", "") == new_sent.replace(" ", "")

                    prot_dict = util.prot_str2entity(old_sent, new_sent, str_prot_set)
                    # 生成输入， 带蛋白质的输入
                    augmented_input = util.gen_input_seq(new_sent, prot_dict)

                    sent_json_form["text"] = augmented_input
                    sent_json_form["event"] = "<extra_id_0> <extra_id_1>"

                    writer.write(json.dumps(sent_json_form) + "\n")
                    sent_json_form.clear()

                    old_sent = ""
                    str_prot_set.clear()

                elif line.startswith("#"):
                    str_prot_set.append(line)
                    prot_info = line.strip().split()

                    assert sent_idx == int(prot_info[3][1:])
                elif line.startswith("$") or line.startswith("@") or line.startswith("%"):
                    continue
                else:
                    old_sent = line
                    sent_idx += 1

if __name__ == '__main__':
    abs_path = "/home/fangsu/myworks/EventGeneration/"
    train_data = "data/raw_data/2013_data/2013_train/"
    train_genia = "data/raw_data/2013_genia_parse/2013_train/"
    train_processed = "data/text2tree/2013_data/2013_train/" # 转换为增广后的数据
    span_train_processed = "data/text2tree/2013_data_span/2013_train/"
    read_data(abs_path+train_data, abs_path+train_genia, abs_path+train_processed, abs_path+span_train_processed)

    devel_data = "data/raw_data/2013_data/2013_devel/"
    devel_genia = "data/raw_data/2013_genia_parse/2013_devel/"
    devel_processed = "data/text2tree/2013_data/2013_devel/"
    span_devel_processed = "data/text2tree/2013_data_span/2013_devel/"
    read_data(abs_path+devel_data, abs_path+devel_genia, abs_path+devel_processed, abs_path+span_devel_processed)

    test_data = "data/raw_data/2013_data/2013_test/"
    test_genia = "data/raw_data/2013_genia_parse/2013_test/"
    test_processed = "data/text2tree/2013_data/2013_test/"
    read_test_data(abs_path+test_data, abs_path+test_genia, abs_path+test_processed)
