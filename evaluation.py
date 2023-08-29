import argparse
import os
from extraction.event_schema import EventSchema
from extraction.predict_parser.tree_predict_parser import TreePredictParser
import process_data.data_utils as util
import utils.eval_utils as eval_utils
import utils.locate_position as locate_position

def read_origin_context(test_file):
    with open(test_file, 'r') as txt_reader:
        context = txt_reader.read()
    return context

def read_file(file_name):
    return [line.strip() for line in open(file_name).readlines()]

def get_trig_idx(origin_a1):
    with open(origin_a1, "r") as reader:
        prot_lines = reader.readlines()
        prot_num = len(prot_lines)
        if prot_num == 0:
            return 0
    return prot_num+1

'''
将test_data中预测的event整理成标准格式(可以上传到评测平台的形式)
'''
# T18	Negative_regulation 0 15	Down-regulation
def trig_normalized_format(sent_trig_list):
    str_trig_dict = dict()
    for trig_entity in sent_trig_list:
        str_trig = trig_entity.trig_idx + "\t" + trig_entity.trig_type + " " + str(trig_entity.context_start) + " " + str(trig_entity.context_end) + "\t" + trig_entity.trig_name
        str_trig_dict[trig_entity.trig_idx] = str_trig

    return str_trig_dict

# 在生成事件时， 尽量补全first_argu_idx, second_argu_idx
# E10	Gene_expression:T27 Theme:T11
# 删除重复的事件， 不能在这里删除， 因为event_idx不一样，所以这里重复事件删除不了
# 另外还涉及到了 event_idx 连续的问题，所以这里不能删除重复事件
def event_normalized_format(sent_event_dict):
    str_event_dict = dict()
    for event_idx, event_entity in sent_event_dict.items():
        trig_info = event_entity.event_type + ":" + event_entity.event_trig_idx
        assert event_entity.first_argu_idx != ""
        first_argu_info = event_entity.first_argu_type + ":" + event_entity.first_argu_idx

        if event_entity.second_argu != None:
            if event_entity.second_argu_type == "Theme":
                second_argu_info = "Theme2" + ":" + event_entity.second_argu_idx
            elif event_entity.second_argu_type == "Cause":
                second_argu_info = event_entity.second_argu_type + ":" + event_entity.second_argu_idx

            str_event = event_idx + "\t" + trig_info + " " + first_argu_info + " " + second_argu_info
            str_event_dict[event_idx] = str_event
        else:
            str_event = event_idx + "\t" + trig_info + " " + first_argu_info
            str_event_dict[event_idx] = str_event

    return str_event_dict


# 返回句子 和 prot_dict
# 用了 yield
def generate_sentence(origin_filename, genia_parse_filename):
    with open(origin_filename, "r") as data_fp, open(genia_parse_filename, "r") as parse_fp:
        old_sent = ""
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
                yield old_sent, new_sent, prot_dict

                old_sent = ""
                str_prot_set.clear()

            elif line.startswith("#"):
                str_prot_set.append(line)

            elif line.startswith("$") or line.startswith("@") or line.startswith("%"):
                continue
            else:
                old_sent = line

def match_sublist(the_list, to_match):
    """
    :param the_list: [1, 2, 3, 4, 5, 6, 1, 2, 4, 5]
    :param to_match: [1, 2]
    :return:
        [(0, 1), (6, 7)]
    """
    len_to_match = len(to_match)
    matched_list = list()
    for index in range(len(the_list) - len_to_match + 1):
        if to_match == the_list[index:index + len_to_match]:
            matched_list += [(index, index + len_to_match - 1)]
    return matched_list

'''
针对一句话进行处理
pred_event_instance['pred_record']: 一句话可能有多个事件
E2 Positive_regulation:T31 Theme:E3 Cause:T6
# record 是一个事件
 考虑的问题：
    (1) 删除重复的trigger
    (2) 删除重复的event
'''
# 补全 trigger的trig_idx, 补全event的event_idx
# 思考重复事件如何解决，涉及到 event_idx 的赋值 问题
def record_to_offset2(test_name, old_sent, new_sent, prot_dict, instance, trig_digit, event_digit, context):
    print(test_name)
    print(old_sent)

    # 这个文件里的重复事件不删除
    # if test_name == "PMC-1359074-03-Results-02.txt" and old_sent.startswith("To investigate whether the interaction with CtBP2 is required"):
    #     print("####################################")
    if old_sent.startswith("The changes in the composition of IL-2Rs were accompanied by inhibition of IL-2-induced"):
        print("###########################")
    if len(prot_dict) == 0:
        return {}, {}, trig_digit, event_digit

    sent_trig_list = list()  # 当前句子中所有的trigger
    sent_event_dict = dict() # 当前句子中所有的event

    # -- 获取trig_name的集合 --
    candidate_trig_list = []
    for record in instance['pred_record']:
        trig_name, trig_matched_list = eval_utils.record2trig(new_sent, record)
        if len(trig_matched_list) > 0:
            candidate_trig_list.append(trig_name)
    #--- 获取 prot_name 的集合
    protname_list = eval_utils.get_protname_set(prot_dict)

    # 生成的事件去重， 生成事件排序，将简单事件放在前面，其实是binding事件， 最后是复杂事件
    retain_event_records = eval_utils.abandon_duplicate_event(instance, protname_list, candidate_trig_list)
    # 生成事件，应该以protein为底
    for record in retain_event_records: # instance['pred_record'], 一个句子里的事件集合, record 是指的一个事件
        trig_entity, event_entity, event_digit = eval_utils.record2event(old_sent, new_sent, prot_dict, record, event_digit, protname_list)
        # 要判断 trigger 是否已经生成，是否已经在sent_trig_dict中存在
        flag = eval_utils.trig_isexist(sent_trig_list, trig_entity)
        if flag == False:
            sent_trig_list.append(trig_entity)

        if event_entity != None:
            assert event_entity.event_idx != ""
            sent_event_dict[event_entity.event_idx] = event_entity

    # 处理 sent_trig_list, 为每个trigger.trig_idx赋值
    trig_dict, trig_digit = eval_utils.assignment_trig_idx(sent_trig_list, trig_digit)
    eval_utils.trig_context_posit(test_name, old_sent, new_sent, prot_dict, trig_dict, context)

    # sent_event_dict中每个event的event_trig_idx赋值
    # 有问题
    eval_utils.complete_event_info(sent_event_dict, trig_dict)

    # 部分事件的论元是Trigger， 将其替换成Trigger引发的事件
    sent_event_dict = eval_utils.event_argu_replace(trig_dict, sent_event_dict)

    # 替换， 嵌套事件用idx 替换
    str_trig_dict = trig_normalized_format(sent_trig_list)
    str_event_dict = event_normalized_format(sent_event_dict)

    return str_trig_dict, str_event_dict, trig_digit, event_digit

# 对 预测的test转变成标准形式，上传到平台测试
# python evaluation.py -g <data-folder-path> -r <offset-folder-path> -p <model-folder-path> -f <data-format>
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-r', dest='offset_folder', default="/home/fangsu/Downloads/bioNlp2011/2011_test_data/") # 官网给出， protein给定了位置位置
    parser.add_argument('-p', dest='pred_folder', default="data/text2tree/2011_pred_test_result/")  # 模型预测出来的事件集合（句子级）
    parser.add_argument('-o', dest='origin_folder', default="data/raw_data/2011_data_one/2011_test/") # offset_file带有offset的文件
    parser.add_argument('-g', dest='genia_parse_folder', default="data/raw_data/2011_genia_parse_one_split_bio/2011_test/") # 分词后的位置
    parser.add_argument('-w', dest='output_folder', default="data/text2tree/2011_standard_format/")  # 转换 2011_pred_test_result 中的文件作为标准格式

    args = parser.parse_args()

    label_schema = EventSchema.read_from_file(
        filename=os.path.join('data/text2tree/2011_data', 'event.schema')
    )

    pred_reader = TreePredictParser(label_constraint=label_schema)
    test_files = os.listdir(args.origin_folder)
    for test_filename in test_files:
        origin_context = read_origin_context(args.offset_folder+test_filename) # 作为最终的检验
        gold_offset_file = args.offset_folder + os.path.splitext(test_filename)[0]+".a1"  # 不带offset 的gold标签，字符序列
        origin_file = args.origin_folder + test_filename
        genia_parse_file = args.genia_parse_folder + test_filename
        pred_file = args.pred_folder + test_filename

        trig_digit = get_trig_idx(gold_offset_file)

        # Read gold event annotation with offsets.
        # 这里需要测试 [(old_sent, new_sent, prots),], 感觉需要 zip 下
        # old_sent_list, new_sent_list, prots_list = [event for event in generate_sentence(origin_file, genia_parse_file)]
        doc_sent_info = [event for event in generate_sentence(origin_file, genia_parse_file)]
        old_sent_list, new_sent_list, prots_list = zip(*doc_sent_info)
        # formed_list是list, 里面的元素是dict
        formed_list, _ = pred_reader.decode(
            gold_list=[],
            pred_list=read_file(pred_file),
            text_list=new_sent_list,
        )

        doc_trig_dict = dict()  # 全局的， 针对一个篇章, 目的是写出
        doc_event_dict = dict()  # 全局的
        event_digit = 1
        for old_sent, new_sent, prot_dict, event_instance in zip(old_sent_list, new_sent_list, prots_list, formed_list): #,针对当前句子
            sent_trig_dict, sent_event_dict, trig_digit, event_digit = record_to_offset2(test_filename, old_sent, new_sent, prot_dict, event_instance, trig_digit, event_digit, origin_context)
            doc_trig_dict.update(sent_trig_dict)
            doc_event_dict.update(sent_event_dict)

        # 将doc_trig_list写出，将doc_event_list写出
        output_filename = os.path.splitext(test_filename)[0] + ".a2"
        with open(args.output_folder + output_filename, "w") as writer:
            for _, trig_seq in doc_trig_dict.items():
                writer.write(trig_seq + "\n")

            for _, str_event in doc_event_dict.items():
                writer.write(str_event + "\n")

if __name__ == "__main__":
    main()