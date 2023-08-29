from process_data.entity import Protein, Trigger, Event
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

    simp_event_set = list()
    bind_event_set = list()
    regu_event_set = list()
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

        evet_type = evet.event_type
        if evet_type in SIMP_TYPE:
            simp_event_set.append(evet)
        elif evet_type in BIND_TYPE:
            bind_event_set.append(evet)
        elif evet_type in REGU_TYPE:
            regu_event_set.append(evet)

    return event_dict

'''
将 string 形式的 trigger 信息转换为实体
@ T18 Negative_regulation S1 0 15 0 15 Down-regulation
# (trig_idx, trig_name, trig_type, trig_oldchar_start, trig_oldchar_end, trig_newchar_start, trig_newchar_end, trig_start, trig_end):
trig_dict中key是trigger名, value是trig_entity
该方法不用于处理多个词组成的trigger
'''
def process_strtrig_one(old_sent, new_sent, str_trig_set):
    trig_dict = dict()
    multi_trigid_set =[]
    for str_trig in str_trig_set:
        trig_info = str_trig.split()
        specialTrigType = ["Entity", 'Anaphora', 'Ubiquitination', 'Protein_modification']
        if len(trig_info) <= 9:
            if trig_info[2] not in specialTrigType:
                new_start, new_end, word_startId, word_endId = relocate(old_sent, new_sent, int(trig_info[4]), int(trig_info[5]))
                trig_entity = Trigger(trig_info[1], trig_info[8], trig_info[2], int(trig_info[4]), int(trig_info[5]), new_start, new_end, word_startId, word_endId)
                trig_dict[trig_info[1]] = trig_entity

        elif len(trig_info) > 9:
            if trig_info[2] not in specialTrigType:
                multi_trigid_set.append(trig_info[1])
    return trig_dict, multi_trigid_set

'''
一个词的trigger 和 多次词的trigger都要进行处理
'''
specialTrigType = ["Entity", 'Anaphora', 'Ubiquitination', 'Protein_modification']
def process_strtrig_two(old_sent, new_sent, str_trig_set):
    trig_dict = dict()
    for str_trig in str_trig_set:
        trig_info = str_trig.split()
        if trig_info[2] not in specialTrigType:
            new_start, new_end, word_startId, word_endId = relocate(old_sent, new_sent, int(trig_info[4]), int(trig_info[5]))
            trig_entity = Trigger(trig_info[1], ' '.join(trig_info[8:]), trig_info[2], int(trig_info[4]), int(trig_info[5]), new_start, new_end, word_startId, word_endId)
            trig_dict[trig_info[1]] = trig_entity
    assert len(trig_dict) == len(str_trig_set)
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

'''
生成序列：Protein和trigger类型以及NoTrigger
prot_trig
'''
def gen_seq(new_sent, prot_dict, trig_dict):
    sent_list = new_sent.split()
    sent_len = len(sent_list)
    seq_result = ["NoTP"]*sent_len # "NoTP"

    prott_len = 0
    trigg_len = 0
    for protein in prot_dict.values():
        prot_start = protein.prot_start #单词的开始位置
        prot_end = protein.prot_end

        prott_len += prot_end-prot_start+1

        for i in range(prot_start, prot_end+1):
            seq_result[i] = "Prot"

        assert sent_len == len(seq_result)

    for trigger in trig_dict.values():
        trig_start = trigger.trig_start
        trig_end = trigger.trig_end
        trig_type = trigger.trig_type
        trigg_len += (trig_end-trig_start)+1
        #assert trig_name == t_name
        for j in range(trig_start, trig_end+1):
            if seq_result[j] == "Prot":
                type = seq_result[j]
                seq_result[j] = type+"_Trig"  #应该是Prot_Trig
            else:
                seq_result[j] = trig_type

        assert sent_len == len(sent_list)

    return seq_result

#--------------------------------------------------------------------------
# 生成输入
def gen_input_seq(sentence, prot_dict):
    prot_entity_dict = {}
    for prot_idx, prot_entity in prot_dict.items():
        prot_entity_dict[prot_idx] = (prot_entity.prot_start, prot_entity.prot_end, "Protein", prot_entity.prot_name)

    sorted_prot_entities = sorted(prot_entity_dict.items(), key=lambda x: x[1][0], reverse=True)
    augmented_input = sentence
    for prot_entity in sorted_prot_entities:
        prot_start, prot_end, prot_type, prot_name = prot_entity[1][0], prot_entity[1][1], prot_entity[1][2], prot_entity[1][3]
        augmented_sent_list = augmented_input.split()
        assert prot_name == " ".join(augmented_sent_list[prot_start:prot_end + 1])
        prefix_str = augmented_sent_list[0:prot_start]
        suffix_str = augmented_sent_list[prot_end + 1:]
        entity_modify = ["<extra_id_0>", prot_name, "<extra_id_2>", prot_type, "<extra_id_1>"]
        augmented_input = " ".join(prefix_str + entity_modify + suffix_str)
    assert len(prot_dict) == augmented_input.split().count("<extra_id_0>")

    return augmented_input


def gen_output_seq(file_name, sent, trig_dict, event_dict):
    if len(trig_dict) == 0 and len(event_dict) == 0:
        return "<extra_id_0> <extra_id_1>"
    elif len(trig_dict) > 0 and len(event_dict) == 0:
        trig_record_set = []
        for trig_idx, trig_entity in trig_dict.items():
            str_trig_info = "<extra_id_0> " + trig_entity.trig_type + " " + trig_entity.trig_name + " <extra_id_1>"
            trig_record_set.append(str_trig_info)

        return " ".join(["<extra_id_0>"] + trig_record_set + ["<extra_id_1>"])

    else:
        tree_event_set = [] # key 是 trigger 的 start, value是event_representation
        for event_idx, event_entity in event_dict.items():
            # trigger information
            event_trig = event_entity.event_trig
            trig_start = event_trig.trig_start
            event_type = event_entity.event_type
            trig_name = event_trig.trig_name
            trig_information = [event_type, trig_name]
            # role information
            first_argu = event_entity.first_argu
            first_argu_type = event_entity.first_argu_type
            first_argu_name = ""
            if isinstance(first_argu, Protein):
                first_argu_name = first_argu.prot_name
            elif isinstance(first_argu, Event):
                argu_trig = first_argu.event_trig
                first_argu_name = argu_trig.trig_name
            else:
                print("process_data/data_utils.py中第344行")
            first_role = ["<extra_id_0>", first_argu_type, first_argu_name, "<extra_id_1>"]

            second_argu = event_entity.second_argu
            second_argu_type = event_entity.second_argu_type
            second_role = []
            if second_argu != None:
                second_argu_name = ""
                if isinstance(second_argu, Protein):
                    second_argu_name = second_argu.prot_name
                elif isinstance(second_argu, Event):
                    argu_trig = second_argu.event_trig
                    second_argu_name= argu_trig.trig_name
                else:
                    print("data_utils.py中第352行")

                second_role = ["<extra_id_0>", second_argu_type, second_argu_name, "<extra_id_1>"]

            event_repre = " ".join(["<extra_id_0>"] + trig_information + first_role + second_role + ["<extra_id_1>"])
            tree_event_set.append((trig_start, event_repre))

        # event_dict 按照 key 值进行排序
        tree_event_set.sort(key=lambda x: (x[0], x[1]))
        result = [value for key, value in tree_event_set]
        assert result.count("<extra_id_0>") == result.count("<extra_id_1>")
        #result1 = " ".join(result).split()
        #print(result1.count("<extra_id_0>"), result1.count("<extra_id_1>"))

    return " ".join(["<extra_id_0>"] + result +["<extra_id_1>"])


def gen_output_seq_span(file_name, sent, trig_dict, event_dict):
    if len(trig_dict) == 0 and len(event_dict) == 0:
        return "<extra_id_0> <extra_id_1>"

    if len(trig_dict) > 0 and len(event_dict) == 0:
        trig_record_set = []
        for trig_idx, trig_entity in trig_dict.items():
            str_trig_info = "<extra_id_0> " + trig_entity.trig_type + " " + trig_entity.trig_name + " <extra_id_1>"
            trig_record_set.append(str_trig_info)

        return " ".join(["<extra_id_0>"] + trig_record_set + ["<extra_id_1>"])

    tree_event_set = []
    for event_idx, event_entity in event_dict.items():
        # trigger information
        event_trig = event_entity.event_trig
        trig_start = event_trig.trig_start
        event_type = event_entity.event_type
        trig_name = event_trig.trig_name
        trig_information = ["<extra_id_0>", event_type, trig_name, "<extra_id_1>"]
        # role information
        first_argu = event_entity.first_argu
        first_argu_type = event_entity.first_argu_type
        first_argu_name = ""
        if isinstance(first_argu, Protein):
            first_argu_name = first_argu.prot_name
        elif isinstance(first_argu, Event):
            argu_trig = first_argu.event_trig
            first_argu_name = argu_trig.trig_name
        else:
            print("process_data/data_utils.py中第344行")
        first_role = ["<extra_id_0>", first_argu_type, first_argu_name, "<extra_id_1>"]

        second_argu = event_entity.second_argu
        second_argu_type = event_entity.second_argu_type
        second_role = []
        if second_argu != None:
            second_argu_name = ""
            if isinstance(second_argu, Protein):
                second_argu_name = second_argu.prot_name
            elif isinstance(second_argu, Event):
                argu_trig = second_argu.event_trig
                second_argu_name= argu_trig.trig_name
            else:
                print("data_utils.py中第352行")

            second_role = ["<extra_id_0>", second_argu_type, second_argu_name, "<extra_id_1>"]

        event_repre = " ".join(trig_information + first_role + second_role)
        tree_event_set.append((trig_start, event_repre))

    # event_dict 按照 key 值进行排序
    tree_event_set.sort(key=lambda x: (x[0], x[1]))
    result = [value for key, value in tree_event_set]
    assert result.count("<extra_id_0>") == result.count("<extra_id_1>")

    return " ".join(["<extra_id_0>"] + result + ["<extra_id_1>"])

##########################################################
# 位置相同， trigger类型不同
# @ T18 Negative_regulation S1 0 15 0 15 Down-regulation
def compute_repeat_trigger(str_trig_set):
    trig_locate_type = {}
    for str_trig in str_trig_set:
        trig_info = str_trig.split()
        trig_start, trig_end, trig_type = trig_info[4], trig_info[5], trig_info[2]

        if trig_start+"_"+trig_end not in trig_locate_type.keys():
            trig_locate_type[trig_start+"_"+trig_end] = trig_type
        else:
            other_trig_type = trig_locate_type[trig_start+"_"+trig_end]
            assert trig_type != other_trig_type
            return True
    return False

# 同一句中 trigger词相同， 但位置不同
# 只考虑trig_name, 不考虑trig_type是否相同
def compute_trigname_same_num(str_trig_set):
    trig_name_set = []
    for str_trig in str_trig_set:
        trig_info = str_trig.split()
        trig_name = " ".join(trig_info[8:])

        if trig_name not in trig_name_set:
            trig_name_set.append(trig_name)
        else:
            return True
    return False

# (3) 同一句话中， protein名称相同，位置不同的句子个数：
# T1 Protein S1 0 5 0 5 BMP-6
def compute_protein_num(str_prot_set):
    protname_list = list()
    for str_prot in str_prot_set:
        prot_info = str_prot.split()
        prot_name = " ".join(prot_info[8:])
        if prot_name not in protname_list:
            protname_list.append(prot_name)
        else:
            return True

    return False

# (3) 同一个序列词， 有些是蛋白质， 有些不被标记为蛋白质
def compute_uncertain_protein(prot_entity_dict, sent):
    prot_name_locate = dict()
    for prot_idx, prot_entity in prot_entity_dict.items():
        prot_start, prot_end = prot_entity.prot_start, prot_entity.prot_end
        prot_name = prot_entity.prot_name
        if prot_name not in prot_name_locate.keys():
            prot_name_locate[prot_name] = [(prot_start, prot_end)]
        else:
            prot_name_locate[prot_name] += [(prot_start, prot_end)]

    token_list = sent.split()
    for prot_name, prot_locate_list in prot_name_locate.items():
        prot_name_list = prot_name.split()
        matched_list = match_sublist(token_list, prot_name_list)
        for matched in matched_list:
            if matched not in prot_locate_list:
                return True
    return False


# (3) 同一个序列词， 有些是trigger， 有些不被标记为trigger
def compute_uncertain_trigger(sent, trigger_entity_dict):
    trig_name_locate = dict()
    for trig_idx, trig_entity in trigger_entity_dict.items():
        trig_start, trig_end = trig_entity.trig_end, trig_entity.trig_end
        trig_name = trig_entity.trig_name
        if trig_name not in trig_name_locate.keys():
            trig_name_locate[trig_name] = [(trig_start, trig_end)]
        else:
            trig_name_locate[trig_name] += [(trig_start, trig_end)]

    token_list = sent.split()
    for trig_name, trig_locate_list in trig_name_locate.items():
        trig_name_list = trig_name.split()
        matched_list = match_sublist(token_list, trig_name_list)
        for matched in matched_list:
            if matched not in trig_locate_list:
                return True
    return False


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







