import copy

from entity import Protein, Trigger, Event
'''
将 string 形式的 protein 信息转换为实体
# T1 Protein S1 19 49 19 49 interferon regulatory factor 4
(prot_idx, prot_name, prot_oldchar_start, prot_oldchar_end, prot_newchar_start, prot_newchar_end, prot_start, prot_end):
'''
def prot_str2entity(old_sent, new_sent, str_prot_set):
    prot_dict = dict()
    for str_prot in str_prot_set:
        prot_info = str_prot.split()
        new_start, new_end, word_start, word_end = relocate(old_sent, new_sent, int(prot_info[4]), int(prot_info[5]))
        prot_entity = Protein(prot_info[1], " ".join(prot_info[8:]), int(prot_info[4]), int(prot_info[5]), new_start, new_end, word_start, word_end)
        prot_entity.prot_context_start = int(prot_info[6])
        prot_entity.prot_context_end = int(prot_info[7])
        prot_entity.prot_mark = prot_info[3]
        prot_dict[prot_info[1]] = prot_entity
    return prot_dict

'''
@ T20 Regulation S2 160 172 306 318 deregulation
(trig_idx, trig_name, trig_type, trig_oldchar_start, trig_oldchar_end, trig_newchar_start, trig_newchar_end, trig_start, trig_end)
'''
def trig_str2entity(old_sent, new_sent, str_trig_set):
    trig_dict = dict()
    for str_trig in str_trig_set:
        trig_info = str_trig.split()
        new_start, new_end, word_start, word_end = relocate(old_sent, new_sent, int(trig_info[4]), int(trig_info[5]))
        trig_entity = Trigger(trig_info[1], " ".join(trig_info[8:]), trig_info[2], int(trig_info[4]), int(trig_info[5]), new_start, new_end, word_start, word_end)
        trig_entity.trig_context_start = int(trig_info[6])
        trig_entity.trig_context_end = int(trig_info[7])
        trig_entity.trig_mark = trig_info[3]
        trig_dict[trig_info[1]] = trig_entity
    return trig_dict

# 根据trigger或protein在old_sent中的位置确定在new_sent中的位置，并且确定出在new_sent中的单词位置
def relocate(old_sent, new_sent, old_start, old_end):
    # 字符的位置
    new_start = locating(old_sent, new_sent, old_start, is_start=True)
    new_end = locating(old_sent, new_sent, old_end, is_start=False)
    assert old_sent[old_start:old_end].replace(" ", "") == new_sent[new_start:new_end].replace(" ", "")

    new_list = new_sent.split()
    # 单词的位置
    if new_start != 0 and new_sent[new_start-1] != " ":
        word_startId = len(new_sent[:new_start].split())-1
    else:
        word_startId = len(new_sent[:new_start].split())

    word_endId = len(new_sent[:new_end].split())-1
    # if new_sent[new_start-1] != " " and new_sent[new_end] != " ":
    #     assert " ".join(new_list[word_startId:word_endId + 1])[:new_end - new_start] == new_sent[new_start:new_end]
    #
    # elif new_sent[new_end] != " ":
    #     assert " ".join(new_list[word_startId:word_endId+1])[:new_end-new_start] == new_sent[new_start:new_end]
    # elif new_start != 0 and new_sent[new_start-1] != " ":
    #     print(new_end-new_start)
    #     print(" ".join(new_list[word_startId:word_endId + 1])[-(new_end - new_start):])
    #     assert " ".join(new_list[word_startId:word_endId + 1])[-(new_end - new_start):] == new_sent[new_start:new_end]
    assert word_startId <=word_endId
    return new_start, new_end, word_startId, word_endId

# 根据字符在old_sent中的位置确定在new_sent中的位置
# is_start是标示目前是计算单词的开始字符位置还是结束字符位置
def locating(old_sent, new_sent, locate, is_start):
    old_str = old_sent[:locate]
    new_str = new_sent[:locate]
    old_str1 = old_sent[:locate].replace(" ", "")
    new_str1 = new_sent[:locate].replace(" ", "")
    new_locate = locate
    if old_str1 == new_str1 and old_str != new_str:
        new_locate += 1

    while old_str1 != new_str1:
        new_locate += 1
        new_str1 = new_sent[:new_locate].replace(" ", "")

    new_str = new_sent[:new_locate]
    if old_str.endswith(" ") and not new_str.endswith(" "):
        new_locate += 1
        new_str1 = new_sent[:new_locate].replace(" ", "")

    if is_start == True:
        if new_sent[new_locate:new_locate+1] == " " and new_sent[:new_locate].replace(" ", "") == new_sent[:new_locate+1].replace(" ", ""):
            new_locate += 1
            new_str1 = new_sent[:new_locate].replace(" ", "")

    assert old_str1 == new_str1
    return new_locate

'''
对source event作处理，处理后可能存在重复，这里对重复事件不作处理
此处没有将Theme2修改称Theme
有一个小小的问题, 有些binding事件可能带不止两个论元
% E11 Positive_regulation:T60 Theme:E9 Cause:T19
'''
def process_strevent(str_event): #注意后面是否有"\n"
    event_info = str_event.strip().split(' ')
    newevent_idx = event_info[1]
    newevent = ''
    if len(event_info) == 4: #事件只有一个论元
        newevent = str_event

    elif len(event_info) == 5:
        if event_info[4].startswith("Theme") or event_info[4].startswith("Cause"):
            newevent = str_event

        else:
            newevent = event_info[0] + " " + event_info[1] + " " + event_info[2] + " " + event_info[3]
            newevent1 = " ".join(event_info[0:4])
            assert newevent == newevent1

    elif len(event_info) >= 6: # 考虑到Binding事件可能不止两个论元,所以
        newevent = event_info[0] + " " + event_info[1] + " " + event_info[2]
        for argu_info in event_info[3:]:
            argu_type = argu_info.split(":")[0]
            argu_idx = argu_info.split(":")[1]
            if argu_type.startswith("Cause"):
                newevent += " " + argu_type+":"+argu_idx
            elif argu_type.startswith("Theme"):
                newevent += " Theme:"+argu_idx
    return newevent, newevent_idx

## % E14 Gene_expression:T31 Theme:T14

SIMP_TYPE = ['Gene_expression', 'Transcription', 'Phosphorylation', 'Protein_catabolism', 'Localization']
BIND_TYPE = ['Binding']
REGU_TYPE = ['Regulation', 'Positive_regulation', 'Negative_regulation']

def event_str2entity(strevent_set, trig_dict, prot_dict):
    event_dict = dict()
    for evet in strevent_set:
        evet_info = evet.split()
        etrig_info = evet_info[2].split(":")
        fargu_info = evet_info[3].split(":")
        if len(evet_info) == 4:
            event_entity = Event()
            event_entity.add_basic_info(evet_info[1], etrig_info[1], etrig_info[0], fargu_info[1], fargu_info[0]) #(event_idx, event_trig_idx, event_type, first_argu_idx, first_argu_type)
            event_dict[evet_info[1]] = event_entity

        elif len(evet_info) == 5:
            event_entity = Event()
            event_entity.add_basic_info(evet_info[1], etrig_info[1], etrig_info[0], fargu_info[1], fargu_info[0])
            sargu_info = evet_info[4].split(":")
            if sargu_info[0].startswith("Theme"):
                event_entity.second_argu_idx = sargu_info[1]
                event_entity.second_argu_type = "Theme"

            elif sargu_info[0].startswith("Cause"):
                event_entity.second_argu_idx = sargu_info[1]
                event_entity.second_argu_type = sargu_info[0]
            else:
                print(sargu_info[0], "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")
            event_dict[evet_info[1]] = event_entity

        else: # 有一个情况, binding带了多个论元(不止两个)
            event_entity = Event()
            ## 第一个论元赋值
            event_entity.add_basic_info(evet_info[1], etrig_info[1], etrig_info[0], fargu_info[1], fargu_info[0])

            ## 第二个论元赋值
            sargu_info = evet_info[4].split(":")
            if sargu_info[0].startswith("Theme"):
                event_entity.second_argu_idx = sargu_info[1]
                event_entity.second_argu_type = "Theme"

            elif sargu_info[0].startswith("Cause"):
                event_entity.second_argu_idx = sargu_info[1]
                event_entity.second_argu_type = sargu_info[0]
            else:
                print(sargu_info[0], "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$")

            ## 多于两个论元的其他论元
            for other_argu in evet_info[5:]:
                other_arguinfo = other_argu.split(":")
                other_argu_idx = other_arguinfo[1]
                other_argu_type = other_arguinfo[0]
                event_entity.other_argu_info[other_argu_idx] = other_argu_type ## 怎么记录论元呢

            event_dict[evet_info[1]] = event_entity

    # 添加event的trig信息以及论元信息
    for event_idx, evet in event_dict.items():
        etrig_idx = evet.event_trig_idx
        etrig = trig_dict[etrig_idx]
        evet.event_trig = etrig
        fargu_idx = evet.first_argu_idx
        if fargu_idx.startswith('T'):
            fargu = prot_dict[fargu_idx]
            evet.first_argu = fargu
        elif fargu_idx.startswith('E'):
            fargu = event_dict[fargu_idx]
            evet.first_argu = fargu

        sargu_idx = evet.second_argu_idx
        if sargu_idx != '' and sargu_idx.startswith('T'):
            sargu = prot_dict[sargu_idx]
            evet.second_argu = sargu
        elif sargu_idx != '' and sargu_idx.startswith('E'):
            sargu = event_dict[sargu_idx]
            evet.second_argu = sargu

        ## other_argu_entity
        if len(evet.other_argu_info) > 0:
            for argu_idx, argu_type in evet.other_argu_info.items(): ## 论元一定是蛋白质
                oargu = prot_dict[argu_idx]
                evet.other_argu_entity[argu_idx] = oargu

    return event_dict

'''
该方法针对的是process_data_styleone.py
生成实体(protein和trigger)和关系一体的增广文本
（1）先生成实体的文本，然后在把关系插入进去
    start_posit = {'S:'+trig_type:(trig_start, trig_end), '1O:'+argu1_type: (argu1_start, argu1_end), '2O:'+argu2_type: (argu2_start, argu2_end)}
    sorted_result = sorted(start_posit.items(), key=lambda x: x[1][0], reverse=True)
'''
# % E16 Negative_regulation:T33 Theme:E17
# def gen_output_text(sentence, prot_dict, trig_dict, event_dict):
#     other_prot_dict = copy.deepcopy(prot_dict)
#     other_trig_dict = copy.deepcopy(trig_dict)
#     for event_idx, event in event_dict.items():
#         first_argu_idx = event.first_argu_idx
#         first_argu_type = event.first_argu_type
#         if first_argu_idx.startswith("T"):
#             prot_entity =
#
#         elif first_argu_idx.startswith("E"):


    # entity_dict = {}
    # for prot_idx, prot_entity in prot_dict.items():
    #     entity_dict[prot_idx]= (prot_entity.prot_start, prot_entity.prot_end, "Protein", prot_entity.prot_name, prot_idx)
    #
    # for trig_idx, trig_entity in trig_dict.items():
    #     entity_dict[trig_idx] = (trig_entity.trig_start, trig_entity.trig_end, trig_entity.trig_type, trig_entity.trig_name, trig_idx)
    #
    # sorted_entities = sorted(entity_dict.items(), key=lambda x: x[1][0], reverse=True)
    #
    # for idx, entity_info in sorted_entities.items():
    #     entity_start, entity_end, entity_type, entity_name = entity_info[0],entity_info[1], entity_info[2], entity_info[3]
    #     assert entity_name == sentence[entity_start:entity_end]
    #     prefix_str = sentence[0:entity_start]
    #     suffix_str = sentence[entity_end:]
    #     entity_modify = [entity_name | entity_type]
    #     sentence = " ".join([prefix_str, entity_modify, suffix_str])

'''
该方法是process_data_styletwo.py中再使用
example:
Therefore, we investigated whether IRF-4 promoter methylation or mutation may be involved in the regulation of IRF-4 expression in leukemia cells.
# T6 Protein S4 35 40 659 664 IRF-4
# T7 Protein S4 111 116 735 740 IRF-4
$ T7 Protein S4 111 116 735 740 IRF-4
@ T21 Regulation S4 97 107 721 731 regulation
@ T22 Gene_expression S4 117 127 741 751 expression
% E4 Regulation:T21 Theme:E5
% E5 Gene_expression:T22 Theme:T7
------------------------------------------------------
trigger识别的输入： ... whether [ IRF-4 | Protein ] promoter ... regulation of [ IRF-4 | Protein ] expression in leukemia cells .
trigger识别的输出： ... in the [ regulation | Regulation ] of IRF-4 [ expression | Gene_expression ] in leukemia cells.
---------------------
关系识别的输入： ... in the [ regulation | Regulation ] of IRF-4 expression in leukemia cells.
关系识别的输出： ... in the regulation of IRF-4 [ expression | Gene_expression | Theme ] in leukemia cells.
处理的是train和devel数据， 如果句子没有蛋白质，忽略该句子
针对T5：
"[" : <extra_id_0> 
"|" : <extra_id_2>
"]" : <extra_id_1>
'''
def gen_augmented_sentence(sentence, prot_dict, trig_dict, event_dict):
    sent_augmented_info = dict()
    # 生成 trigger的输入
    prot_entity_dict = {}
    for prot_idx, prot_entity in prot_dict.items():
        prot_entity_dict[prot_idx] = (prot_entity.prot_start, prot_entity.prot_end, "Protein", prot_entity.prot_name, prot_idx)

    sorted_prot_entities = sorted(prot_entity_dict.items(), key=lambda x: x[1][0], reverse=True)
    augmented_trig_input = sentence
    for prot_entity in sorted_prot_entities:
        prot_start, prot_end, prot_type, prot_name = prot_entity[1][0], prot_entity[1][1], prot_entity[1][2], prot_entity[1][3]
        augmented_sent_list = augmented_trig_input.split()
        assert prot_name == " ".join(augmented_sent_list[prot_start:prot_end+1])
        prefix_str = augmented_sent_list[0:prot_start]
        suffix_str = augmented_sent_list[prot_end+1:]
        entity_modify = ["<extra_id_0>", prot_name, "<extra_id_2>",  prot_type, "<extra_id_1>"]
        augmented_trig_input = " ".join(prefix_str + entity_modify + suffix_str)

    # 生成 trigger 的增广输出
    trig_entity_dict = {}
    for trig_idx, trig_entity in trig_dict.items():
        trig_entity_dict[trig_idx] = (trig_entity.trig_start, trig_entity.trig_end, trig_entity.trig_type, trig_entity.trig_name, trig_idx)

    sorted_trig_entities = sorted(trig_entity_dict.items(), key=lambda x: x[1][0], reverse=True)

    augmented_trig_output = sentence
    for trig_entity in sorted_trig_entities:
        trig_start, trig_end, trig_type, trig_name = trig_entity[1][0],trig_entity[1][1], trig_entity[1][2], trig_entity[1][3]
        augmented_trig_list = augmented_trig_output.split(" ")
        assert trig_name == " ".join(augmented_trig_list[trig_start:trig_end+1])
        prefix_str = augmented_trig_list[0:trig_start]
        suffix_str = augmented_trig_list[trig_end+1:]
        entity_modify =["<extra_id_0>", trig_name, "<extra_id_2>", trig_type, "<extra_id_1>"]
        augmented_trig_output = " ".join(prefix_str+entity_modify+suffix_str)
    sent_augmented_info["trigger"] = (augmented_trig_input, augmented_trig_output)
    #--------------------------------------------------------------
    # 针对每个trigger，生成对应的语料
    # relation的输入：
    reverse_sorted_trig_entities = sorted(trig_entity_dict.items(), key=lambda x: x[1][0], reverse=False) # 升序
    augmented_relation_set = []
    for trig_entity in reverse_sorted_trig_entities:
        trig_idx = trig_entity[0]
        trig_start = trig_entity[1][0]
        trig_end = trig_entity[1][1]
        trig_type = trig_entity[1][2]
        trig_name = trig_entity[1][3]
        augmented_relation_input = sentence
        augmented_relation_inputlist = augmented_relation_input.split()
        prefix_str = augmented_relation_inputlist[0:trig_start]
        suffix_str = augmented_relation_inputlist[trig_end+1:]
        trig_modify = ["<extra_id_0>", trig_name, "<extra_id_2>", trig_type, "<extra_id_1>"]
        augmented_relation_input = " ".join(prefix_str + trig_modify + suffix_str)

        # 考虑一个trigger引发多个事件的情况
        related_arguments = {}
        for event_idx, event in event_dict.items():
            if trig_idx == event.event_trig_idx:
                first_argu_idx = event.first_argu_idx
                first_argu = event.first_argu
                if first_argu_idx.startswith("T"):
                    related_arguments[first_argu_idx] = (first_argu.prot_start, first_argu.prot_end, event.first_argu_type, first_argu)
                elif first_argu_idx.startswith("E"):
                    argu_event_trig = event_dict[first_argu_idx].event_trig
                    related_arguments[first_argu_idx] = (argu_event_trig.trig_start, argu_event_trig.trig_end, event.first_argu_type, argu_event_trig)

                second_argu_idx  = event.second_argu_idx
                second_argu = event.second_argu
                if second_argu != None:
                    if second_argu_idx.startswith("T"):
                        related_arguments[second_argu_idx] = (second_argu.prot_start, second_argu.prot_end, event.second_argu_type, second_argu)
                    elif second_argu_idx.startswith("E"):
                        second_argu_event_trig = event_dict[second_argu_idx].event_trig
                        related_arguments[second_argu_idx] = (second_argu_event_trig.trig_start, second_argu_event_trig.trig_end, event.second_argu_type, second_argu_event_trig)

        # 对论元集合排序
        augmented_relation_out = sentence
        sorted_relation_entities = sorted(related_arguments.items(), key=lambda x: x[1][0], reverse=True)  #降序
        for entity_info in sorted_relation_entities: # entity_info 可能 trigger, protein
            entity_start, entity_end, argu_type, entity = entity_info[1][0], entity_info[1][1], entity_info[1][2], entity_info[1][3]
            augmented_relation_outlist = augmented_relation_out.split()
            prefix_str = augmented_relation_outlist[0: entity_start]
            suffix_str = augmented_relation_outlist[entity_end+1:]
            entity_modify = ""
            if isinstance(entity, Protein):
                entity_modify = ["<extra_id_0>", entity.prot_name, "<extra_id_2> Protein <extra_id_2>", argu_type, "<extra_id_1>"]
            elif isinstance(entity, Trigger):
                entity_modify = ["<extra_id_0>", entity.trig_name, "<extra_id_2>", entity.trig_type, "<extra_id_2>", argu_type, "<extra_id_1>"]
            elif isinstance(entity, Event):
                argu_trig = entity.event_trig
                entity_modify = ["<extra_id_0>", argu_trig.trig_name, "<extra_id_2>", argu_trig.trig_type, "<extra_id_2>", argu_type, "<extra_id_1>"]
                print("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
            assert entity_modify != ""
            augmented_relation_out = " ".join(prefix_str + entity_modify + suffix_str)

        # 一个trigger的论元数据结束
        augmented_relation_set.append((augmented_relation_input, augmented_relation_out))
    sent_augmented_info["relation"] = augmented_relation_set
    return sent_augmented_info

'''
针对T5
"[" : <extra_id_0> 
"|" : <extra_id_2>
"]" : <extra_id_1>
'''
def gen_test_augmented_sentence(sentence, prot_dict):
    augmented_sents = {}
    prot_entity_dict = {}
    for prot_idx, prot_entity in prot_dict.items():
        prot_entity_dict[prot_idx]= (prot_entity.prot_start, prot_entity.prot_end, "Protein", prot_entity.prot_name, prot_idx)

    # 如果句子中没有标注的 protein， 则输入就是原句子
    sorted_prot_entities = sorted(prot_entity_dict.items(), key=lambda x: x[1][0], reverse=True)
    augmented_prot_sentence = sentence
    for prot_entity in sorted_prot_entities:
        prot_start, prot_end, prot_type, prot_name = prot_entity[1][0], prot_entity[1][1], prot_entity[1][2], prot_entity[1][3]
        augmented_sent_list = augmented_prot_sentence.split()
        assert prot_name == " ".join(augmented_sent_list[prot_start:prot_end+1])
        prefix_str = augmented_sent_list[0:prot_start]
        suffix_str = augmented_sent_list[prot_end+1:]
        entity_modify = ["<extra_id_0>", prot_name, "<extra_id_2>",  prot_type, "<extra_id_1>"]
        augmented_prot_sentence = " ".join(prefix_str + entity_modify + suffix_str)
    augmented_sents['trigger'] = (augmented_prot_sentence, "")
    augmented_sents['relation'] = []
    return augmented_sents
