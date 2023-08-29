'''
统计 trian 和 devel 所有 句子中：
（1）同一个trigger带有两个不同类型的个数有多少个
（2）同一句中 trigger词相同， 但位置不同
'''
import os
import data_utils as util
import delete_nostandard as dns
def read_data(data_path, genia_parse_path):
    files = os.listdir(data_path)
    total_sent_num = 0
    sent_hastrig_num = 0
    sent_hasevent_num = 0
    sent_hasprot_num = 0
    sent_sametrig_difftrigtype = 0
    sent_trigname_num = 0
    sent_sameprotein_num = 0 # 相同的蛋白质，位置不同
    sent_uncertain_protein_num = 0 # 有些子串在句子中被标记为protein, 有些不被标记为蛋白质
    sent_uncertain_trigger_num = 0 # 子串在句子中有些被标记为trigger, 有些不被标记为trigger


    for file_name in files:
        with open(data_path + file_name, "r") as data_fp, open(genia_parse_path + file_name, "r") as genia_parse_fp:
            old_sent = ""
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
                    assert old_sent.replace(" ", "") == new_sent.replace(" ", "")

                    # -----------------------------(1)删除跨句子事件--------------------------------------------
                    valid_event_set = dns.delete_invalid_event_two(str_event_set, str_prot_idxs, str_event_idx)

                    # 将str_prot, str_trig, str_event转换成event
                    prot_dict = util.prot_str2entity(old_sent, new_sent, str_prot_set)
                    # jprot_dict = util.prot_str2entity(old_sent, new_sent, str_jprot_set)
                    trig_dict = util.trig_str2entity(old_sent, new_sent, str_trig_set)
                    event_dict = util.event_str2entity(valid_event_set, trig_dict, prot_dict)
                    assert len(valid_event_set) == len(event_dict)
                    total_sent_num += 1

                    if len(str_trig_set) > 0:
                        sent_hastrig_num += 1
                    if len(valid_event_set) > 0:
                        sent_hasevent_num += 1
                    if len(str_prot_set) > 0:
                        sent_hasprot_num += 1

                    # （1）同一个trigger带有两个不同类型的个数有多少个
                    result1 = util.compute_repeat_trigger(str_trig_set)
                    if result1:
                        sent_sametrig_difftrigtype += 1

                    # （2） 同名 trigger 的句子个数， 位置不同
                    result2 = util.compute_trigname_same_num(str_trig_set)
                    if result2:
                        sent_trigname_num += 1

                    # 同一个句子中， 一些是trigger 词， 另一些不是trigger词
                    result5 = util.compute_uncertain_trigger(new_sent, trig_dict)
                    if result5:
                        sent_uncertain_trigger_num += 1

                    # (3) protein名称相同，位置不同的句子个数：
                    result3 = util.compute_protein_num(str_prot_set)
                    if result3:
                        sent_sameprotein_num += 1

                    # (3) 同一个序列词， 有些是蛋白质， 有些不被标记为蛋白质
                    result4 = util.compute_uncertain_protein(prot_dict, new_sent)
                    if result4:
                        sent_uncertain_protein_num += 1

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
                    #assert sent_idx == int(prot_info[3][1:])

                elif line.startswith("$"):
                    str_jprot_set.append(line)

                elif line.startswith("@"):# @ T20 Regulation S2 160 172 306 318 deregulation
                    trig_info = line.split()
                    if trig_info[2] not in ["Entity", 'Anaphora', 'Ubiquitination', 'Protein_modification']:
                        str_trig_set.append(line)

                elif line.startswith("%"):
                    str_event, event_idx = util.process_strevent(line)
                    if str_event not in str_event_set: ## 删除重复event
                        str_event_set.append(str_event)
                        str_event_idx.append(event_idx)

                else:
                    old_sent = line

    print("句子总数：", total_sent_num)
    print("有trigger的句子个数：", sent_hastrig_num)
    print("有event的句子个数：", sent_hasevent_num)
    print("相同位置的trigger, 不同的类型的句子个数：", sent_sametrig_difftrigtype)
    print("相同的trigger_name的句子个数：", sent_trigname_num)
    print("标记了蛋白质的句子个数：", sent_hasprot_num)
    print("相同的蛋白质名不同的位置句子个数", sent_sameprotein_num)
    print("相同的字符串一些被标记为蛋白质，一些不被标记为蛋白质的句子个数", sent_uncertain_protein_num)
    print("相同的字符串一些被标记为trigger，一些不被标记为trigger的句子个数", sent_uncertain_trigger_num)
    print("#########################################################")


if __name__ == '__main__':
    abs_path = "/home/fangsu/myworks/EventGeneration/"
    train_data = "data/raw_data/2011_data_one/2011_train/"
    train_genia = "data/raw_data/2011_genia_parse_one_split_bio/2011_train/"
    read_data(abs_path + train_data, abs_path + train_genia)
    print("#########################################################")

    devel_data = "data/raw_data/2011_data_one/2011_devel/"
    devel_genia = "data/raw_data/2011_genia_parse_one_split_bio/2011_devel/"
    read_data(abs_path + devel_data, abs_path + devel_genia)
    print("%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%")