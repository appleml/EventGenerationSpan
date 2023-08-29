import os
from process_data import delete_nostandard as dns
import data_utils as data_utils
import json

#instance与#
def data_event_check(origin_data, parse_data, processed_data):
    instance_num = 0
    instance_event_num = 0
    sent_set = []
    with open(processed_data, "r") as json_fp:
        for line in json_fp:
            instance_num += 1
            if line.strip() != "":
                instance = json.loads(line.strip())

                if instance["event"] == "<extra_id_0>  <extra_id_1>":
                    instance_event_num +=1
                sent_set.append(instance["text"])


    total_sent_num = 0
    no_event_one_num = 0
    no_event_two_num = 0
    event_sent_num = 0
    valid_event_num = 0
    number = 0 # 句子中事件为0而trigger不为零的个数
    files = os.listdir(origin_data)
    sent_number = 0
    for file_name in files:
        with open(origin_data + "/" + file_name, "r") as data_fp, open(parse_data +"/"+file_name, "r") as genia_parse_fp:
            old_sent = ""
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

                    genia_parse_line = genia_parse_fp.readline()
                    genia_parse_line = genia_parse_line.strip().strip('\n')
                    while len(genia_parse_line) != 0:
                        genia_parse = genia_parse_line.split("\t###\t")
                        # genia_line信息处理，为识别trigger做准备
                        # (word+"\t"+pos+"\t"+fea1+"\t"+fea2+"\t"+prot+"\t"+cjprot+"\t"+ctrig+"\t"+jprot+"\t"+trigtype+"\n")
                        genia_line = genia_parse[0]
                        genia_info = genia_line.split("\t")
                        gword_set.append(genia_info[0])

                        genia_parse_line = genia_parse_fp.readline()
                        genia_parse_line = genia_parse_line.strip().strip('\n')

                    new_sent = " ".join(gword_set)
                    assert old_sent.replace(" ", "") == new_sent.replace(" ", "")

                    # --------(1)删除掉四个类型的trigger(2)删除重复事件(3)删除跨句子事件(3)删除多词trigger引发的事件--------------
                    valid_event_set = dns.delete_invalid_event_two(str_event_set, str_prot_idxs, str_event_idx)
                    # if len(valid_event_set) == 0:
                    #     no_event_one_num += 1
                    # -----（4）一个trigger词有多个类型， 保留一个类型(保留简单类型)
                    valid_trig_set, valid_event_set = dns.delete_invalid_event(str_trig_set, valid_event_set)
                    valid_event_num += len(valid_event_set)
                    # ----------------------------------------------------------
                    if new_sent not in sent_set:
                        print(file_name)
                        print(new_sent)
                        print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@2")
                        sent_number += 1
                    # if len(str_event_set) == 0 and len(valid_event_set) == 0:
                    #     no_event_one_num += 1
                    #
                    # elif len(str_event_set) != 0 and len(valid_event_set) == 0:
                    #     no_event_two_num += 1
                    # else:
                    #     assert len(str_event_set) != 0 and len(valid_event_set) != 0
                    #     event_sent_num += 1
                    #
                    # if len(valid_event_set) == 0 and (len(str_trig_set) != 0 or len(str_prot_set) != 0):
                    #     number += 1

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

                elif line.startswith("@"):  # @ T20 Regulation S2 160 172 306 318 deregulation
                    trig_info = line.split()
                    if trig_info[2] not in ["Entity", 'Anaphora', 'Ubiquitination', 'Protein_modification']:
                        str_trig_set.append(line)

                elif line.startswith("%"):
                    str_event, event_idx = data_utils.process_strevent(line)
                    if str_event not in str_event_set:  ## 删除重复event
                        str_event_set.append(str_event)
                        str_event_idx.append(event_idx)
                else:
                    old_sent = line
                    total_sent_num += 1

    print("转换JSON文件后， 样例个数：",instance_num,)
    #print("转为JSON文件后，事件不为空的样例个数：", instance_event_num)
    # print("不包含事件的句子个数：", no_event_one_num + no_event_two_num)
    # print("句子中有prot 或者 trigger而没有事件的句子个数：", number)
    # print("有效的事件个数：", valid_event_num)
    print("语料中分割的句子个数：", total_sent_num)
    # print("包含事件的句子的个数：", event_sent_num)
    # print("句子中原本就没有标注事件", no_event_one_num)
    # print("句子有标注的事件，不合格删除了", no_event_two_num)
    print(sent_number)


def compute_label_event_number(path):
    total_event_number = 0
    files = os.listdir(path)
    for file_name in files:
        filename, extension = os.path.splitext(file_name)
        if extension == ".a2":
            with open(path+"/"+file_name, "r") as origin_fp:
                for line in origin_fp:
                    if line.startswith("E"):
                        total_event_number += 1

    print("语料中标注的事件总数：", total_event_number)

if __name__ == "__main__":
    absolute_path = "//"
    # label_data = "/home/fangsu/Downloads/bioNlp2011/2011_train_data"
    # compute_label_event_number(label_data)

    origin_train_data = absolute_path + "data/2011_data_one/2011_train"
    parse_train_data = absolute_path + "data/2011_genia_parse_one_split_bio/2011_train"
    processed_train_data = absolute_path + "data/2011_tree_subtype/train.json"
    data_event_check(origin_train_data, parse_train_data, processed_train_data)

    print("##################################################")

    # label_devel_data = "/home/fangsu/Downloads/bioNlp2011/2011_devel_data"
    # compute_label_event_number(label_devel_data)

    # origin_devel_data = absolute_path + "data/2011_data_one/2011_devel"
    # parse_devel_data = absolute_path + "data/2011_genia_parse_one_split_bio/2011_devel"
    # processed_devel_data = absolute_path + "data/2011_tree_subtype/devel.json"
    # data_event_check(origin_devel_data, parse_devel_data, processed_devel_data)