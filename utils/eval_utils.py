import copy
from process_data.entity import Protein, Trigger, Event

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
返回 (1) trig_name
    (2) len(trig_matched_list) > 0
'''
def record2trig(new_sent, record):
    token_list = new_sent.split()
    trig_name = record['trigger']
    trig_matched_list = match_sublist(token_list, trig_name.split())
    return trig_name, trig_matched_list

# 保证一个句子中总有一个事件的落脚点是protein
# 假设 origin_trigname_list中是重复的
# 思路：一次从origin_trigname_list中删除一个
def remove_duplicate_event(retain_event_record, protname_list, origin_trigname_list):
    # 不和规矩的事件就删除掉
    wellformed_events = list()
    for event_record in retain_event_record:
        event_type = event_record['type']
        # 如果论元为空，则可以生成 event
        if len(event_record['roles']) != 0:
            flag = False
            for _, argu_type, argu_name in event_record['roles']:
                if argu_type == "Cause":
                    flag = True
                    break
            if event_type in ['Gene_expression', 'Transcription', 'Phosphorylation', 'Protein_catabolism', 'Localization', 'Binding']:
                if flag == False:
                    wellformed_events.append(event_record)

            elif event_type in ['Regulation', 'Positive_regulation', 'Negative_regulation']:
                wellformed_events.append(event_record)

        else:
            wellformed_events.append(event_record)

    #------------------------------------------------------------------------------------------------
    prev_event_num = len(wellformed_events)
    epoch = 0
    while len(wellformed_events) > 0 and (epoch == 0 or prev_event_num != len(wellformed_events)):
        epoch += 1
        prev_event_num = len(wellformed_events)
        event_record_list = list()
        for event_record in wellformed_events:
            trig_name = event_record['trigger']
            flag = False
            for _, argu_type, argu_name in event_record['roles']:
                if argu_name not in protname_list+origin_trigname_list:
                    flag = True
                    break
            if flag == False:
                event_record_list.append(event_record)
            else:
                origin_trigname_list.remove(trig_name)

        wellformed_events = copy.deepcopy(event_record_list)

    return wellformed_events


# 生成的事件去重，有问题，如句子， "(ii) Transfection and expression of EBNA-LP alone had no effect on LMP1 expression "
# 将简单事件放在前面，其实是binding事件， 最后是复杂事件
# 确保论元不是protein,就是trigger
'''
(1) 去除重复出现的事件
(2) simple事件和binding事件的论元必须是protein， 否则删除
(3) simple 事件的类型必须是Theme， 不能是Cause
(4) complex事件必须一个论元是Theme， 另一个是Cause
'''
def abandon_duplicate_event(instance, protname_list, trigname_list):
    # 删除重复的事件
    retain_event_record = list()
    str_event_info = list()
    for event_record in instance['pred_record']:
        event_type = event_record['type']
        trig_name = event_record['trigger']
        argu_info = ""
        for _, argu_type, argu_name in event_record['roles']:
            argu_info += argu_type + "_" + argu_name+"_"
        if argu_info.endswith("_"):
            argu_info = argu_info[0:-1]
        str_event = event_type+"_"+trig_name+"_"+argu_info
        if str_event not in str_event_info:
            str_event_info.append(str_event)
            retain_event_record.append(event_record)
    # 删除不合适的事件， 比如 Simple 事件有Cause 论元

    # 某些事件的论元可能既不是protein, 也不是trig, 则必须删除
    retain_event_record = remove_duplicate_event(retain_event_record, protname_list, trigname_list)

    # 如果 retain_event_record中 所有的事件没有以蛋白质为根基，则有问题
    flag = False
    for event_record in retain_event_record:
        prot_flag = list()
        for _, argu_type, argu_name in event_record['roles']:
            if argu_name in protname_list:
                prot_flag.append(True)
            else:
                prot_flag.append(False)
        if len(list(set(prot_flag))) == 1 and prot_flag[0] == True: # 说明有以蛋白质为论元的事件
            flag = True
            break
    if flag == False:
        return []

    # 事件排序, 另外要保证simple, binding的论元是protein
    new_trigname_list = copy.deepcopy(trigname_list)
    simple_record = list()
    binding_record = list()
    complex_prot_record = list() # 如果仅有一个论元， 论元是protein， 如果有两个论元，都是蛋白质才可以
    complex_trig_record = list()
    for event_record in retain_event_record:
        event_type = event_record['type']
        trig_name = event_record['trigger']
        # 如果论元为空，则可以生成 event
        if len(event_record['roles']) == 0:
            if event_type in ['Gene_expression', 'Transcription', 'Phosphorylation', 'Protein_catabolism', 'Localization']:
               simple_record.append(event_record)
            elif event_type in ['Binding']:
                binding_record.append(event_record)
            elif event_type in ['Regulation', 'Positive_regulation', 'Negative_regulation']:
               complex_prot_record.append(event_record)

        else: # 论元不为空， 再次检查论元以及论元类型，合适就存储
            argu_name_flag = []
            argu_name_list = []
            argu_type_list = []
            flag = True
            for _, argu_type, argu_name in event_record['roles']:
                argu_name_list.append(argu_name)
                argu_type_list.append(argu_type)
                if argu_name not in protname_list + new_trigname_list:
                    flag = False
                    new_trigname_list.remove(trig_name)
                    break
                if argu_name in protname_list:
                    argu_name_flag.append(True)
                else:
                    assert argu_name in new_trigname_list
                    argu_name_flag.append(False)

            if flag == True:
                if event_type in ['Gene_expression', 'Transcription', 'Phosphorylation', 'Protein_catabolism', 'Localization']:
                    # 要确保论元类型也是Theme， 论元名字是protein
                    p_flag = True
                    for argu_name in argu_name_list:
                        if argu_name not in protname_list:
                            p_flag = False
                            break
                    if p_flag == True and "Cause" not in argu_type_list:
                        simple_record.append(event_record)
                    else:
                        new_trigname_list.remove(trig_name)

                elif event_type in ['Binding']:
                    p_flag = True
                    for argu_name in argu_name_list:
                        if argu_name not in protname_list:
                            p_flag = False
                            break
                    if p_flag == True and argu_type_list[0] == "Theme":
                        binding_record.append(event_record)

                elif event_type in ['Regulation', 'Positive_regulation', 'Negative_regulation']:
                    if len(argu_type_list) == 2 and "Cause" not in argu_type_list:
                        break
                    if len(list(set(argu_name_flag))) == 1 and argu_name_flag[0] == True:
                        complex_prot_record.append(event_record)
                    else:
                        if trig_name not in argu_name_list:  # 如果trigger和论元trigger是同一个词，
                            complex_trig_record.append(event_record)


    return simple_record+binding_record+complex_prot_record+complex_trig_record

'''
record 转化为一个 event 实体
argu_list: (argu_type, argu_name, argu_matched_list)
return: 
    (1) Trigger
    (2) None or Event
    (3) trig_idx, event_idx
在该函数中， 对生成Trigger不编码（即不设定trig_idx, 一个trigger可能引发两个事件，如果对每个事件的trigger都编码 trig_idx, 可能同一个trigger词有两个trig_idx）
# 现在不进行trig_idx, 因为会有重复， 多个record可能是同一个trigger
如果论元是protein, 就赋值first_argu_idx, second_argu_idx
如果论元不是protein而是trigger， 那么该trigger必须出现在candidate_trig_list中，否则record就废除不考虑其能构建事件
'''
def record2event(old_sent, new_sent, prot_dict, record, event_digit, protname_list):
    #probably_arguname_list = protname_list + trigname_list
    token_list = new_sent.split()
    event_type = record['type']
    trig_name = record['trigger']
    trig_matched_list = match_sublist(token_list, trig_name.split())
    if len(trig_matched_list) == 0:
        return None, None, event_digit

    argu_list = []  # 论元可能是一个或者两个
    for _, argu_type, argu_name in record['roles']:
        argu_matched_list = match_sublist(token_list, argu_name.split())

        if len(argu_matched_list) == 0:
             sys.stderr.write("[Cannot reconstruct]: %s %s\n" % (argu_name, token_list))
        else:
            if (argu_type, argu_name, argu_matched_list) not in argu_list:
                argu_list += [(argu_type, argu_name, argu_matched_list)]

    event = None
    if len(trig_matched_list) == 1:
        trig_start, trig_end = trig_matched_list[0]
        oldchar_start, oldchar_end, newchar_start, newchar_end = complete_position(old_sent, new_sent, trig_name, trig_start, trig_end)
        trig_entity = Trigger("", trig_name, event_type, oldchar_start, oldchar_end, newchar_start, newchar_end, trig_start, trig_end)
        if len(argu_list) == 0:
            return trig_entity, None, event_digit

        elif len(argu_list) == 1:  # 该事件可能有一个论元
            argu_type, argu_name, argu_matched_list = argu_list[0]
            # if argu_name not in probably_arguname_list or argu_name == trig_name:
            #     return trig_entity, None, event_digit

            assert argu_type == "Theme"
            if len(argu_matched_list) == 1:
                argu_entity = get_argu_entity(old_sent, new_sent, prot_dict, argu_name, argu_matched_list[0][0], argu_matched_list[0][1])
                if argu_entity != None:
                    event, event_digit = gen_event(event_digit, event_type, trig_entity, argu_entity, argu_type)
                    return trig_entity, event, event_digit

                return trig_entity, None, event_digit

            elif len(argu_matched_list) > 1:  # 找一个离trigger最近的论元
                argu_entity = get_best_argu(old_sent, new_sent, prot_dict, trig_start, trig_end, argu_name, argu_matched_list)
                if argu_entity != None:
                    event, event_digit = gen_event(event_digit, event_type, trig_entity, argu_entity, argu_type)
                    return trig_entity, event, event_digit
                return trig_entity, None, event_digit

        elif len(argu_list) == 2:  # trigger可能有两个论元
            # argu_list += [(argu_type, argu_name, argu_matched_list)]
            first_argu_type, first_argu_name, first_matched_list = argu_list[0]
            second_argu_type, second_argu_name, second_matched_list = argu_list[1]

            # if first_argu_name not in probably_arguname_list or second_argu_name not in probably_arguname_list:
            #     return trig_entity, None, event_digit

            if first_argu_name == trig_name or second_argu_name == trig_name:
                return trig_entity, None, event_digit

            if len(first_matched_list) == 1 and len(second_matched_list) == 1:
                first_argu_start, first_argu_end = first_matched_list[0]
                first_argu_entity = get_argu_entity(old_sent, new_sent, prot_dict, first_argu_name, first_argu_start, first_argu_end)
                second_argu_start, second_argu_end = second_matched_list[0]
                second_argu_entity = get_argu_entity(old_sent, new_sent, prot_dict, second_argu_name, second_argu_start, second_argu_end)
                if first_argu_entity != None and second_argu_entity != None:
                    event, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu_entity, first_argu_type, second_argu_entity, second_argu_type)
                    return trig_entity, event, event_digit

                return trig_entity, None, event_digit

            elif len(first_matched_list) == 1 and len(second_matched_list) > 1:
                first_argu_start, first_argu_end = first_matched_list[0]
                first_argu_entity = get_argu_entity(old_sent, new_sent, prot_dict, first_argu_name, first_argu_start, first_argu_end)
                second_argu_entity = get_best_argu(old_sent, new_sent, prot_dict, trig_start, trig_end, second_argu_name, second_matched_list)
                if first_argu_entity != None and second_argu_entity != None:
                    event, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu_entity, first_argu_type, second_argu_entity, second_argu_type)
                    return trig_entity, event, event_digit

                return trig_entity, None, event_digit

            elif len(first_matched_list) > 1 and len(second_matched_list) == 1:
                first_argu_entity = get_best_argu(old_sent, new_sent, prot_dict, trig_start, trig_end, first_argu_name, first_matched_list)
                second_argu_start, second_argu_end = second_matched_list[0]
                second_argu_entity = get_argu_entity(old_sent, new_sent, prot_dict, second_argu_name, second_argu_start, second_argu_end)
                if first_argu_entity != None and second_argu_entity != None:
                    event, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu_entity, first_argu_type, second_argu_entity, second_argu_type)
                    return trig_entity, event, event_digit
                return trig_entity, None, event_digit

            elif len(first_matched_list) > 1 and len(second_matched_list) > 1:
                first_argu_entity = get_best_argu(old_sent, new_sent, prot_dict, trig_start, trig_end, first_argu_name, first_matched_list)
                second_argu_entity = get_best_argu(old_sent, new_sent, prot_dict, trig_start, trig_end, second_argu_name, second_matched_list)
                if first_argu_entity != None and second_argu_entity != None:
                    event, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu_entity, first_argu_type, second_argu_entity, second_argu_type)
                    return trig_entity, event, event_digit
                return trig_entity, None, event_digit
        else:
            return trig_entity, None, event_digit

    elif len(trig_matched_list) > 1:
        if len(argu_list) == 0:
            # 怎么确定trig_entity, 直接取第一个
            trig_start, trig_end = trig_matched_list[0]
            oldchar_start, oldchar_end, newchar_start, newchar_end = complete_position(old_sent, new_sent, trig_name, trig_start, trig_end)
            trig_entity = Trigger("", trig_name, event_type, oldchar_start, oldchar_end, newchar_start, newchar_end, trig_start, trig_end)
            return trig_entity, None, event_digit

        #event_type, trig_name, trig_matched_list
        elif len(argu_list) == 1: # 只有一个论元
            argu_type, argu_name, argu_matched_list = argu_list[0]

            # if argu_name not in probably_arguname_list or argu_name == trig_name:
            #     return None, None, event_digit

            assert argu_type == "Theme"
            if len(argu_matched_list) == 0:
                # 在已有的trig_list中拿到最后一个trigger, 然后从trig_matched_list中一个位置大于 trig_end的
                # trig_matched = get_proper_trig(exist_trig_list, trig_matched_list)
                # trig_entity = gen_trigger(old_sent, new_sent, event_type, trig_name, trig_matched[0], trig_matched[1])
                # return trig_entity, None
                return None, None, event_digit

            elif len(argu_matched_list) == 1:
                trig_entity = get_trig_entity(old_sent, new_sent, event_type, trig_name, trig_matched_list, argu_matched_list[0][0], argu_matched_list[0][1])
                argu_entity = get_argu_entity(old_sent, new_sent, prot_dict, argu_name, argu_matched_list[0][0], argu_matched_list[0][1])

                # 确定一个离protein 最近的trigger词
                event, event_digit = gen_event(event_digit, event_type, trig_entity, argu_entity, argu_type)
                return trig_entity, event, event_digit

            elif len(argu_matched_list) > 1: # 找一个离trigger最近的论元
                trig_entity, event, event_digit = get_trig_argu_entity(old_sent, new_sent, prot_dict, event_digit, event_type, trig_name, trig_matched_list, argu_type, argu_name, argu_matched_list)
                return trig_entity, event, event_digit

        elif len(argu_list) == 2:
            first_argu_type, first_argu_name, first_matched_list = argu_list[0]
            second_argu_type, second_argu_name, second_matched_list = argu_list[1]

            # if first_argu_name not in probably_arguname_list and second_argu_name not in probably_arguname_list:
            #     return None, None, event_digit

            if first_argu_name == trig_name or second_argu_name == trig_name:
                return None, None, event_digit

            if len(first_matched_list) == 1 and len(second_matched_list) == 1:
                trig_entity, event, event_digit = gen_event_multrig_twoargu_11(old_sent, new_sent, prot_dict,
                                                                               event_digit, event_type, trig_name, trig_matched_list,
                                                                               first_argu_name, first_argu_type, first_matched_list,
                                                                               second_argu_name, second_argu_type, second_matched_list)
                return trig_entity, event, event_digit

            elif len(first_matched_list) == 1 and len(second_matched_list) > 1:
                trig_entity, event, event_digit = gen_event_multrig_twoargu_12(old_sent, new_sent, prot_dict,
                                                                               event_digit, event_type, trig_name, trig_matched_list,
                                                                               first_argu_name, first_argu_type, first_matched_list,
                                                                               second_argu_name, second_argu_type, second_matched_list)
                return trig_entity, event, event_digit

            elif len(first_matched_list) > 1 and len(second_matched_list) == 1:
                trig_entity, event, event_digit = gen_event_multrig_twoargu_21(old_sent, new_sent, prot_dict,
                                                                               event_digit, event_type, trig_name, trig_matched_list,
                                                                               first_argu_name, first_argu_type, first_matched_list,
                                                                               second_argu_name, second_argu_type, second_matched_list)

            elif len(first_matched_list) > 1 and len(second_matched_list) > 1:
                trig_entity, event, event_digit = gen_event_multrig_twoargu_22(old_sent, new_sent, prot_dict,
                                                                               event_digit, event_type, trig_name, trig_matched_list,
                                                                               first_argu_name, first_argu_type, first_matched_list,
                                                                               second_argu_name, second_argu_type, second_matched_list)
        else:
            print("eval_utils.py第152行")


        return trig_entity, event, event_digit

'''
给定triggr，first_argu, second_argu,生成事件
'''
def gen_event(event_digit, event_type, event_trig, first_argu, first_argu_type, *args):
    event = Event()
    event.event_idx = "E"+str(event_digit)
    event.event_type = event_type
    event.event_trig = event_trig
    event.first_argu = first_argu
    event.first_argu_type = first_argu_type
    if len(args) > 0:
        assert len(args) == 2
        event.second_argu = args[0]
        event.second_argu_type = args[1]
    return event, event_digit+1

def get_argu_entity(old_sent, new_sent, prot_dict, argu_name, argu_start, argu_end):
    prot_idx = check_isprot(prot_dict, argu_name, argu_start, argu_end)

    if prot_idx != "":
        return prot_dict[prot_idx]
    else:
        trig_entity = gen_trigger(old_sent, new_sent, "", argu_name, argu_start, argu_end)
        return trig_entity

'''
找到匹配类的论元原则
前提条件 len(argu_matched_list) > 1
(1) trigger 的位置给定
(1) 首先判断其是否是蛋白质， 如果是蛋白质， 一定是在标注的蛋白质中选择
(2) 在所有的matched中， 选择离trigger位置最近
matched_text 对应着matched_list
'''
def get_best_argu(old_sent, new_sent, prot_dict, trig_start, trig_end, argu_name, argu_matched_list):
    cand_prot_dict = get_matched_protein(prot_dict, argu_name, argu_matched_list)

    if len(cand_prot_dict) != 0:
        argu_prot = get_nearest_protein(cand_prot_dict, trig_start, trig_end)  # 离 trigger 最近的 protein
        return argu_prot

    # argu_name 不是 protein
    else:
        dist_info = list()
        for argu_matched in argu_matched_list:
            dist = compute_distance(trig_start, trig_end, argu_matched[0], argu_matched[1])
            dist_info += [(dist, argu_matched)]
        dist_info.sort(key=lambda x: (x[0]))
        _, argu_matched = dist_info[0]

        argu_trig = gen_trigger(old_sent, new_sent, "", argu_name, argu_matched[0], argu_matched[1])
        return argu_trig

# 根据 prot_dict得到 prot_name 的集合
def get_protname(prot_dict):
    protname_list = list()
    for _, prot_entity in prot_dict.items():
        prot_name = prot_entity.prot_name
        if prot_name not in protname_list:
            protname_list.append(prot_name)
    return protname_list

def get_best_trig(old_sent, new_sent, event_type, trig_name, trig_matched_list, argu_start, argu_end):
    dist_info_list = []
    for trig_matched in trig_matched_list:
        cand_trig_start, cand_trig_end = trig_matched[0], trig_matched[1]
        dist = compute_distance(cand_trig_start, cand_trig_end, argu_start, argu_end)
        dist_info_list += [(dist, trig_matched)]
    dist_info_list.sort(key=lambda x: (x[0]))
    _, trig_matched = dist_info_list[0]

    trig_entity = gen_trigger(old_sent, new_sent, event_type, trig_name, trig_matched[0], trig_matched[1])
    return trig_entity

'''
判断论元的matched_text是否是protein， 
如果是protein, 则在prot_dict中找，否则就是trigger(即嵌套事件）
# 适合 len(argu_matched_list) == 1 的情况
'''
def check_isprot(prot_dict, text_str, text_start, text_end):
    for prot_idx, prot_entity in prot_dict.items():
        prot_start, prot_end, prot_name = prot_entity.prot_start, prot_entity.prot_end, prot_entity.prot_name
        if prot_start == text_start and prot_end == text_end and prot_name == text_str:
            return prot_idx
    return ""

'''
仅仅判断argu_name是否可能是protein
'''
def is_protein(prot_dict, matched_text) -> object:
    protname_list = get_protname(prot_dict)
    if matched_text in protname_list:
        return True
    return False

'''
有些相同的字符被标注为protein, 有些没有被标注为protein, 
如果matched_text是 prot_name, 要保证从prot_list中选择一个离trigger最近
前提概要： 在同一句话中，同一个字符串， 有些被标注为protein, 有些没有被标注为蛋白质
目的：将没有标注为protein 的 matched 去掉
'''
def get_matched_protein(prot_dict, matched_text, matched_list):
    protname_list = get_protname(prot_dict)
    cand_prot_dict = dict()
    if matched_text in protname_list:
        for prot_idx, prot_entity in prot_dict.items():
            prot_start, prot_end = prot_entity.prot_start, prot_entity.prot_end
            if (prot_start, prot_end) in matched_list:
                cand_prot_dict[prot_idx] = prot_entity

    return cand_prot_dict

'''
(1) 应该确定了matched_text是protein同名
(2) 找到一个与trigger离的最近的protein
'''
import sys
def get_nearest_protein(prot_dict, trig_start, trig_end):
    distance_list = list()
    for prot_idx, prot_entity in prot_dict.items():
        prot_start, prot_end = prot_entity.prot_start, prot_entity.prot_end
        dist = compute_distance(trig_start, trig_end, prot_start, prot_end)
        distance_list += [(dist, prot_idx)]

    distance_list.sort(key=lambda x: (x[0]))
    _, prot_idx = distance_list[0]

    return prot_dict[prot_idx]

'''
 len(tirgger_matched_list) > 1
 len(argu_match) == 1 论元只有一个， 选距离论元最近的trigger
'''
def get_trig_entity(old_sent, new_sent, trig_type, trig_name, trig_matched_list, argu_start, argu_end):
    distance_list = list()
    for matched in trig_matched_list:
        cand_trig_start, cand_trig_end = matched[0], matched[1]
        distance = compute_distance(cand_trig_start, cand_trig_end, argu_start, argu_end)
        distance_list += [(distance, matched)]

    distance_list.sort(key=lambda x: x[0])

    _, (trig_start, trig_end) = distance_list[0]

    trig_entity = gen_trigger(old_sent, new_sent, trig_type, trig_name, trig_start, trig_end)

    return trig_entity


'''
len(argu_list) == 1 只有一个（Theme）论元
len(trig_matched_list) > 1
len(argu_matched_list) > 1
新生成的事件不能和以往的事件重复(预测是不理会)
event_set: 已经生成的事件
判断 argu_name 是否属于protein 或 trigger， 如果都不输入，则不会生成事件
'''
def get_trig_argu_entity(old_sent, new_sent, prot_dict, event_digit, event_type, trig_name, trig_matched_list, argu_type, argu_name, argu_matched_list):
    cand_prot_list = get_matched_protein(prot_dict, argu_name, argu_matched_list)
    if len(cand_prot_list) > 0:
        distance_info = list()
        for trig_matched in trig_matched_list:
            for prot_idx, prot_entity in cand_prot_list.items():
                prot_start, prot_end = prot_entity.prot_start, prot_entity.prot_end
                dist = compute_distance(trig_matched[0], trig_matched[1], prot_start, prot_end)
                distance_info += [(dist, trig_matched, prot_entity)]
        distance_info.sort(key=lambda x: (x[0]))
        _,  trig_matched, prot_entity = distance_info[0]

        trig_entity = gen_trigger(old_sent, new_sent, event_type, trig_name, trig_matched[0], trig_matched[1])
        event, event_digit = gen_event(event_digit, event_type, trig_entity, prot_entity, argu_type)
        return trig_entity, event, event_digit

    else:
        distance_info = list()
        for trig_matched in trig_matched_list:
            for argu_matched in argu_matched_list:
                dist = compute_distance(trig_matched[0], trig_matched[1], argu_matched[0], argu_matched[1])
                distance_info += [(dist, trig_matched, argu_matched)]

        distance_info.sort(key=lambda x: (x[0]))
        _, trig_matched, prot_matched = distance_info[0]

        trig_entity = gen_trigger(old_sent, new_sent, event_type, trig_name, trig_matched[0], trig_matched[1])
        argu_entity = gen_trigger(old_sent, new_sent, "", argu_name, prot_matched[0], prot_matched[1])
        event, event_digit = gen_event(event_digit, event_type, trig_entity, argu_entity, argu_type)
        return trig_entity, event, event_digit

'''
len(trig_matched_list) > 1:
len(first_matched_list) == 1 and len(second_matched_list) == 1:
一共有两个论元
'''
def gen_event_multrig_twoargu_11(old_sent, new_sent, prot_dict,
                                 event_digit, event_type, trig_name, trig_matched_list,
                                 first_argu_name, first_argu_type, first_matched_list,
                                 second_argu_name, second_argu_type, second_matched_list):

    first_matched = first_matched_list[0]
    first_argu = get_argu_entity(old_sent, new_sent, prot_dict, first_argu_name, first_matched[0], first_matched[1])
    if first_argu != None:
        trig_entity = get_best_trig(old_sent, new_sent, event_type, trig_name, trig_matched_list, first_matched[0], first_matched[1])
        second_matched = second_matched_list[0]
        second_argu = get_argu_entity(old_sent, new_sent, prot_dict, second_argu_name, second_matched[0], second_matched[1])
        if second_argu != None:
            event_entity, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu, first_argu_type, second_argu, second_argu_type)
            return trig_entity, event_entity, event_digit
        return trig_entity, None, event_digit
    return None, None, event_digit

'''
len(triger_matched_list) > 1
len(first_matched_list) == 1 and len(second_matched_list) > 1
想法：
根据一个论元的位置，确定trigger的位置， 然后根据trigger的位置再确定另一个论元的位置
'''
def gen_event_multrig_twoargu_12(old_sent, new_sent, prot_dict,
                                 event_digit, event_type, trig_name, trig_matched_list,
                                 first_argu_name, first_argu_type, first_matched_list,
                                 second_argu_name, second_argu_type, second_matched_list):

    first_matched = first_matched_list[0]
    first_argu = get_argu_entity(old_sent, new_sent, prot_dict, first_argu_name, first_matched[0], first_matched[1])
    if first_argu != None:
        # 根据first_argu， 确定trigger的位置
        trig_entity = get_best_trig(old_sent, new_sent, event_type, trig_name, trig_matched_list, first_matched[0], first_matched[1])
        # 根据trigger的位置，确定另一个论元
        second_argu = get_best_argu(old_sent, new_sent, prot_dict, trig_entity.trig_start, trig_entity.trig_end, second_argu_name, second_matched_list)
        if second_argu != None:
            event, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu, first_argu_type, second_argu, second_argu_type)
            return trig_entity, event, event_digit
        return trig_entity, None, event_digit

    return None, None, event_digit


'''
len(triger_matched_list) > 1
len(first_matched_list) > 1 and len(second_matched_list) == 1
想法：
根据一个论元的位置，确定trigger的位置， 然后根据trigger的位置再确定另一个论元的位置
'''
def gen_event_multrig_twoargu_21(old_sent, new_sent, prot_dict,
                                 event_digit, event_type, trig_name, trig_matched_list,
                                 first_argu_name, first_argu_type, first_matched_list,
                                 second_argu_name, second_argu_type, second_matched_list):
    second_matched = second_matched_list[0]
    second_argu = get_argu_entity(old_sent, new_sent, prot_dict, second_argu_name, second_matched[0], second_matched[1])
    if second_argu != None:
        # 根据second_argu， 确定trigger的位置
        trig_entity = get_best_trig(old_sent, new_sent, event_type, trig_name, trig_matched_list, second_matched[0], second_matched[1])

        # 根据trigger的位置，确定另一个论元
        first_argu = get_best_argu(old_sent, new_sent, prot_dict, trig_entity.trig_start, trig_entity.trig_end, first_argu_name, first_matched_list)
        if first_argu != None:
            event, event_digit = gen_event(event_digit, event_type, trig_entity, first_argu, first_argu_type, second_argu, second_argu_type)
            return trig_entity, event, event_digit
        return trig_entity, None, event_digit

    return None, None, event_digit

'''
len(triger_matched_list) > 1
len(first_matched_list) > 1 and len(second_matched_list) > 1
想法：
根据一个论元的位置，确定trigger的位置， 然后根据trigger的位置再确定另一个论元的位置
'''
def gen_event_multrig_twoargu_22(old_sent, new_sent, prot_dict,
                                 event_digit, event_type, trig_name, trig_matched_list,
                                 first_argu_name, first_argu_type, first_matched_list,
                                 second_argu_name, second_argu_type, second_matched_list):

    first_prot_dict = get_matched_protein(prot_dict, first_argu_name, first_matched_list)
    second_prot_dict = get_matched_protein(prot_dict, second_argu_name, second_matched_list)
    distance_information = []
    for trig_matched in trig_matched_list:
        first_dist_info = []
        if len(first_prot_dict) > 0:
            for first_prot_idx, first_prot in first_prot_dict.items():
                first_prot_start, first_prot_end = first_prot.prot_start, first_prot.prot_end
                first_dist = compute_distance(trig_matched[0], trig_matched[1], first_prot_start, first_prot_end)
                first_dist_info += [(first_dist, (first_prot_start, first_prot_end), first_prot_dict[first_prot_idx])]
        else:
            for first_matched in first_matched_list:
                first_dist = compute_distance(trig_matched[0], trig_matched[1], first_matched[0], first_matched[1])
                first_dist_info += [(first_dist, first_matched, None)]

        second_dist_info = []
        if len(second_prot_dict) > 0:
            for second_prot_idx, second_prot in second_prot_dict.items():
                second_prot_start, second_prot_end = second_prot.prot_start, second_prot.prot_end
                second_dist = compute_distance(trig_matched[0], trig_matched[1], second_prot_start, second_prot_end)
                second_dist_info += [(second_dist, (second_prot_start, second_prot_end), second_prot_dict[second_prot_idx])]

        else:
            for second_matched in second_matched_list:
                second_dist = compute_distance(trig_matched[0], trig_matched[1], second_matched[0], second_matched[1])
                second_dist_info += [(second_dist, second_matched, None)]

        for first_dist in first_dist_info:
            first_dist, first_matched, first_entity = first_dist
            for second_dist in second_dist_info:
                second_dist, second_matched, second_entity = second_dist
                distance_information += [(first_dist+second_dist, trig_matched, first_matched, first_entity, second_matched, second_entity)]

    distance_information.sort(key=lambda x:x[0])
    _, trig_matched, first_matched, first_entity, second_matched, second_entity = distance_information[0]

    trig_entity = gen_trigger(old_sent, new_sent, event_type, trig_name, trig_matched[0], trig_matched[1])

    if first_entity != None and second_entity != None:
        event, event_digit = gen_event(event_digit, event_type, trig_entity, first_entity, first_argu_type, second_entity, second_argu_type)
        return trig_entity, event, event_digit

    elif first_entity != None and second_entity == None:
        second_entity = get_argu_entity(old_sent, new_sent, prot_dict, second_argu_name, second_matched[0], second_matched[1])
        event, event_digit = gen_event(event_digit, event_type, trig_entity, first_entity, first_argu_type, second_entity, second_argu_type)
        return trig_entity, event, event_digit

    elif first_entity == None and second_entity != None:
        first_entity = get_argu_entity(old_sent, new_sent, prot_dict, first_argu_name, first_matched[0], first_matched[1])
        event, event_digit = gen_event(event_digit, event_type, trig_entity, first_entity, first_argu_type, second_entity, second_argu_type)
        return trig_entity, event, event_digit

    elif first_entity == None and second_entity == None:
        first_entity = get_argu_entity(old_sent, new_sent, prot_dict, first_argu_name, first_matched[0], first_matched[1])
        second_entity = get_argu_entity(old_sent, new_sent, prot_dict, second_argu_name, second_matched[0], second_matched[1])
        event, event_digit = gen_event(event_digit, event_type, trig_entity, first_entity, first_argu_type, second_entity, second_argu_type)
        return trig_entity, event, event_digit

    return trig_entity, None, event_digit

'''
len(trig_matched_list) > 1
一共有两个论元， 给定了两个论元的offset, 确定trigger的set
'''
def find_best_trig(old_sent, new_sent, trig_name, trig_type, trig_matched_list, fist_argu_matched, second_argu_matched):
    trig_distance_information = []
    first_argu_start, first_argu_end = fist_argu_matched[0], fist_argu_matched[1]
    second_argu_start, second_argu_end = second_argu_matched[0], second_argu_matched[1]
    for trig_matched in trig_matched_list:
        cand_trig_start, cand_trig_end = trig_matched[0], trig_matched[1]
        dist1 = compute_distance(cand_trig_start, cand_trig_end, first_argu_start, first_argu_end)
        dist2 = compute_distance(cand_trig_start, cand_trig_end, second_argu_start, second_argu_end)
        trig_distance_information += [(dist1+dist2, trig_matched)]

    trig_distance_information.sort(key=lambda x : (x[0]))
    _, trig_matched = trig_distance_information[0]

    trig_entity = gen_trigger(old_sent, new_sent, trig_type, trig_name, trig_matched[0], trig_matched[1])

    return trig_entity

'''
计算trigger 和 argument 之间的距离
'''
def compute_distance(trig_start, trig_end, argu_start, argu_end):
    if trig_start == argu_start and trig_end == argu_end:
        return 0
    elif trig_start > argu_end:
        return trig_start-argu_end
    elif argu_start > trig_end:
        return argu_start-trig_end
    else:
        print("eval_utils.py中第326行代码出问题啦")
        sys.exit()

# 确定 trig_oldchar_start, trig_oldchar_end, trig_newchar_start, trig_newchar_end,
def complete_position(old_sent, new_sent, trig_name, trig_start, trig_end):
    token_list = new_sent.split()
    start_prefix = " ".join(token_list[0:trig_start])
    if start_prefix == "":
        trig_newchar_start = len(start_prefix)
    else:
        trig_newchar_start = len(start_prefix) + 1  ## 加1是空格
    end_prefix = " ".join(token_list[0:trig_end + 1])
    trig_newchar_end = len(end_prefix)
    assert new_sent[trig_newchar_start:trig_newchar_end] == trig_name

    # 根据trig_newchar_start, trig_newchar_end， 然后计算 trig_oldchar_start， trig_oldchar_end
    trig_oldchar_start = comp_oldchar_position(old_sent, new_sent, trig_newchar_start, False)
    trig_oldchar_end = comp_oldchar_position(old_sent, new_sent, trig_newchar_end, True)
    # print(file_name, "--------------------------")
    # print(old_sent[trig_entity.trig_oldchar_start:trig_entity.trig_oldchar_end], "---", trig_entity.trig_name)
    assert old_sent[trig_oldchar_start:trig_oldchar_end].replace(" ", "") == trig_name.replace(" ", "")

    return trig_oldchar_start, trig_oldchar_end, trig_newchar_start, trig_newchar_end

# 根据在new_sent中的位置计算在old_sent中的位置
def comp_oldchar_position(old_sent, new_sent, newchar_position, is_end):
    oldchar_position = newchar_position

    #while len(newchar_prefix) != len(oldchar_prefix):
    while new_sent[0:newchar_position].replace(" ", "") != old_sent[0:oldchar_position].replace(" ",""):
        oldchar_position = oldchar_position-1

    if is_end == True and new_sent[0:newchar_position].replace(" ", "") == old_sent[0:oldchar_position-1].replace(" ",""):
        oldchar_position = oldchar_position-1

    assert old_sent[0:oldchar_position].replace(" ", "") == new_sent[0:newchar_position].replace(" ", "")
    return oldchar_position

def gen_trigger(old_sent, new_sent, event_type, trig_name, trig_start, trig_end):
    oldchar_start, oldchar_end, newchar_start, newchar_end = complete_position(old_sent, new_sent, trig_name, trig_start, trig_end)
    trig_entity = Trigger("", trig_name, event_type, oldchar_start, oldchar_end, newchar_start, newchar_end, trig_start, trig_end)

    return trig_entity

'''
 trigger 设定 tirg_idx， 确定 trig_type 
 根据trigger的位置， 设置 trig_idx
 argu_trig的类型 无需补充， 用已经有类型的event替换
 注意： 一个trigger 可以触发多个事件，不要给同一个trigger 赋值不同的值
'''
def complete_trigidx(event_list, trig_digit, event_digit):
    position = list()
    trig_list = list()
    for event in event_list:
        trig_entity = event.event_trig
        trig_list.append(trig_entity)
        position += [(trig_entity.trig_start, event)]

    trig_digit, trig_entity_list = clear_trig(trig_list, trig_digit)

    position.sort(key=lambda x: x[0])
    event_dict = dict()
    for _, event in position:
        event.event_trig.trig_idx = "T"+str(trig_digit)
        event.event_idx = "E" + str(event_digit)
        event_dict["E" + str(event_digit)] = event
        trig_digit += 1
        event_digit += 1

    return trig_digit, event_digit, trig_list, event_dict

'''
去除重复的trigger
相同的trigger就是：trig_start, trig_end, trig_type相同
y
'''
def clear_trig(trig_list, trig_digit):
    trig_entity_list = list()
    trig_info = list()
    for trig in trig_list:
        if (trig.trig_start, trig.trig_end, trig.trig_type, trig.trig_name) not in trig_info:
            trig_info.append((trig.trig_start, trig.trig_end, trig.trig_type, trig.trig_name))
            trig_entity_list.append(trig)

    trig_info.sort(key=lambda x:x[0])
    for trig in trig_entity_list:
        trig.trig_idx = "T"+str(trig_digit)
        trig_digit += 1

    return trig_digit, trig_entity_list

'''
论元是 trigger, 要改成 Event
Trigger(self, trig_idx, trig_name, trig_type, trig_oldchar_start, trig_oldchar_end, trig_newchar_start, trig_newchar_end, trig_start, trig_end):
思路： 分成两份： 一份的论元是protein，
               另一份的论元是论元是trigger, 重点就是修改这一份中的事件，将其trigger论元转换成event论元
'''
def trigargu2eventargu(sent_event_dict):
    complete_events = dict()
    imperfect_events = dict()
    epoch = 0
    while epoch == 0 or len(imperfect_events) != 0:
        imperfect_again = dict()
        for evet in sent_event_dict:
            evet_idx = evet.event_idx
            first_argu = evet.first_argu
            second_argu = evet.second_argu
            if second_argu == None:
                if isinstance(first_argu, (Protein, Event)):
                    complete_events[evet.event_idx] = evet

                elif isinstance(first_argu, Trigger):
                    e_idx = trigargu2eventidx(complete_events, first_argu, evet_idx)
                    if e_idx != "":
                        arguevent_entity = complete_events[e_idx]
                        evet.first_argu = arguevent_entity
                        complete_events[evet.event_idx] = evet  # evet 此时 已经 由trigger论元替换成event论元了

                    else:
                        imperfect_again[evet.event_idx] = evet

            else:
                if isinstance(first_argu, (Protein, Event)) and isinstance(second_argu, (Protein,Event)):
                    complete_events[evet.event_idx] = evet

                elif isinstance(first_argu, (Protein,Event)) and isinstance(second_argu, Trigger):
                    argu_evnet_idx = trigargu2eventidx(complete_events, second_argu, evet_idx)
                    if argu_evnet_idx != "":
                        arguevent_entity = complete_events[argu_evnet_idx]
                        evet.second_argu = arguevent_entity
                        complete_events[evet.event_idx] = evet  # evet 此时 已经 由trigger论元替换成event论元了

                    else:
                        imperfect_again[evet.event_idx] = evet

                elif isinstance(first_argu, Trigger) and isinstance(second_argu, (Protein,Event)):
                    argu_evnet_idx = trigargu2eventidx(complete_events, first_argu, evet_idx)
                    if argu_evnet_idx != "":
                        arguevent_entity = complete_events[argu_evnet_idx]
                        evet.first_argu = arguevent_entity
                        complete_events[evet.event_idx] = evet  # evet 此时 已经 由trigger论元替换成event论元了

                    else:
                        imperfect_again[evet.event_idx] = evet

                elif isinstance(first_argu, Trigger) and isinstance(second_argu, Trigger):
                    first_argu_evnet_idx = trigargu2eventidx(complete_events, first_argu, evet_idx)
                    if first_argu_evnet_idx != "":
                        arguevent_entity = complete_events[first_argu_evnet_idx]
                        evet.first_argu = arguevent_entity

                    second_argu_evnet_idx = trigargu2eventidx(complete_events, second_argu, evet_idx)
                    if second_argu_evnet_idx != "":
                        arguevent_entity = complete_events[second_argu_evnet_idx]
                        evet.first_argu = arguevent_entity

                    if first_argu_evnet_idx != "" and second_argu_evnet_idx != "":
                        complete_events[evet.event_idx] = evet
                    else:
                        imperfect_again[evet.event_idx] = evet

        imperfect_events = copy.deepcopy(imperfect_again)
    return complete_events

'''
argu_trig 是 imcomplete_eventidx的 trigger 论元， 需要将其转换成事件
'''
import utils.metric_utils as metric_utils
def trigargu2eventidx(event_dict, argu_trig, imcomplete_eventidx):
    for event_idx, event_trig in event_dict.items():
        flag = metric_utils.trig_isequal(event_trig, argu_trig)
        if flag == True and event_idx != imcomplete_eventidx:
            return event_idx
    return ""

def get_protname_set(prot_dict):
    protname_list = list()
    for _, prot_entity in prot_dict.items():
        if prot_entity.prot_name not in protname_list:
            protname_list.append(prot_entity.prot_name)
    return protname_list

#--- 定位 trigger 在 篇章中的位置
def trig_context_posit(file_name, old_sent, new_sent, prot_dict, trig_dict, context):
    new_sent_list = new_sent.split()
    #_, first_prot = next(iter(prot_dict.items()))
    first_prot = prot_dict.get(next(iter(prot_dict)))
    for trig_idx, trig_entity in trig_dict.items():
        trig_start = trig_entity.trig_start
        trig_end = trig_entity.trig_end
        trig_name = trig_entity.trig_name
        ## 先计算 trig_newchar_start, trig_newchar_end
        start_prefix = " ".join(new_sent_list[0:trig_start])
        if start_prefix == "":
            trig_entity.trig_newchar_start = len(start_prefix)
        else:
            trig_entity.trig_newchar_start = len(start_prefix)+1 ## 加1是空格
        end_prefix = " ".join(new_sent_list[0:trig_end+1])
        trig_entity.trig_newchar_end = len(end_prefix)
        assert new_sent[trig_entity.trig_newchar_start:trig_entity.trig_newchar_end] == trig_entity.trig_name

        ## 根据trig_newchar_start, trig_newchar_end， 然后计算 trig_oldchar_start， trig_oldchar_end
        trig_entity.trig_oldchar_start = comp_oldchar_position(old_sent, new_sent, trig_entity.trig_newchar_start, False)
        trig_entity.trig_oldchar_end = comp_oldchar_position(old_sent, new_sent, trig_entity.trig_newchar_end, True)
        #print(file_name, "--------------------------")
        #print(old_sent[trig_entity.trig_oldchar_start:trig_entity.trig_oldchar_end], "---", trig_entity.trig_name)
        assert old_sent[trig_entity.trig_oldchar_start:trig_entity.trig_oldchar_end].replace(" ", "") == trig_name.replace(" ", "")

        ## trig_oldchar_start， trig_oldchar_end，最后计算 trig_context_start, trig_context_end
        trigger_name = old_sent[trig_entity.trig_oldchar_start:trig_entity.trig_oldchar_end]
        trig_entity.context_start, trig_entity.context_end = comp_article_position(file_name, old_sent, context, trigger_name, trig_entity.trig_oldchar_start, first_prot)

def comp_article_position(file_name, old_sent, context, trig_name, trig_oldchar_start, first_prot):
    trig_context_start = -1
    trig_context_end = -1
    prot_context_start = first_prot.prot_context_start
    prot_context_end = first_prot.prot_context_end

    assert prot_context_start != -1 and prot_context_end != -1
    prot_oldchar_start = int(first_prot.prot_oldchar_start)
    prot_oldchar_end = int(first_prot.prot_oldchar_end)

    if prot_context_start == prot_oldchar_start and prot_context_end == prot_oldchar_end:
        trig_context_start = trig_oldchar_start
        trig_context_end = trig_oldchar_start+len(trig_name)

    elif trig_oldchar_start > prot_oldchar_end:
        relative_distance = trig_oldchar_start - prot_oldchar_end
        trig_context_start = prot_context_end + relative_distance
        trig_context_end = trig_context_start+len(trig_name)
        #print('---', context[trig_context_start:trig_context_end], '---', trig_name, '---')

    elif trig_oldchar_start == prot_oldchar_end:
        if old_sent[trig_oldchar_start-1:trig_oldchar_start] != " ":
            trig_context_start = prot_context_end
            trig_context_end = trig_context_start + len(trig_name)
        else:
            trig_context_start = prot_context_end-1
            trig_context_end = trig_context_start+len(trig_name)

    elif trig_oldchar_start == prot_oldchar_start:
        trig_context_start = prot_context_start
        trig_context_end += trig_context_start+len(trig_name)+1

    elif trig_oldchar_start+len(trig_name) < prot_oldchar_start:
        relative_distance = prot_oldchar_start - trig_oldchar_start
        trig_context_start = prot_context_start - relative_distance
        trig_context_end += trig_context_start+len(trig_name)+1

    elif prot_oldchar_start < trig_oldchar_start and trig_oldchar_start+len(trig_name) <= prot_oldchar_end:
        relative_distance = trig_oldchar_start - prot_oldchar_start
        trig_context_start = prot_context_start + relative_distance
        trig_context_end = trig_context_start+len(trig_name)

    elif trig_oldchar_start < prot_oldchar_start:
        relative_distance = prot_oldchar_start - trig_oldchar_start
        trig_context_start = prot_context_start-relative_distance
        trig_context_end = trig_context_start+len(trig_name)

    else:
        print("there is a serious problem!")

    if context[trig_context_start:trig_context_end] != trig_name:
        trig_context_start += 1
        trig_context_end += 1

    if context[trig_context_start:trig_context_end] != trig_name:
        print(file_name)
        print(old_sent)
        print(trig_name)
        print(trig_context_start, '---', trig_context_end)
        print('---', context[trig_context_start:trig_context_end], '---', trig_name, '---')
        print("#################$$$$$$$")

    return trig_context_start, trig_context_end

#----------------------------------------------------------------------------------------
from process_data.entity import Event, Trigger, Protein


#比较两个trigger是否相同
def trig_isequal(trig1, trig2):
    if trig1.trig_name == trig2.trig_name and trig1.trig_start == trig2.trig_start and \
                    trig1.trig_end == trig2.trig_end and trig1.trig_type == trig2.trig_type:
        return True
    return False

## 判断新生成的 trigger 是否已存在
def trig_isexist(trig_list, specific_trig):
    if specific_trig == None:
        return True
    for trig in trig_list:
        flag = trig_isequal(trig, specific_trig)
        if flag == True:
            return True
    return False

#-------------------------------------------------------------
#
def find_trigidx(main_trig, argu_trig):
    if argu_trig.trig_type == "":
        if main_trig.trig_name == argu_trig.trig_name and main_trig.trig_start == argu_trig.trig_start and \
                main_trig.trig_end == argu_trig.trig_end:
            return main_trig.trig_idx
    if main_trig.trig_name == argu_trig.trig_name and main_trig.trig_start == argu_trig.trig_start and \
                    main_trig.trig_end == argu_trig.trig_end and main_trig.trig_type == argu_trig.trig_type:
        return main_trig.trig_idx
    return ""

def confirm_trig_idx(trig_dict, argu_trig):
    for trig_idx, trig in trig_dict.items():
        cand_trig_idx = find_trigidx(trig, argu_trig)
        if cand_trig_idx != "":
            return cand_trig_idx
    return ""
# ----------------------------------------------------------------
def assignment_trig_idx(trig_list, trig_digit):
    trig_dict = dict()
    for trig_entity in trig_list:
        trig_entity.trig_idx = "T"+str(trig_digit)
        trig_dict["T"+str(trig_digit)] = trig_entity
        trig_digit += 1
    return trig_dict, trig_digit
'''
# sent_event_dict中的每个event_trig赋值 event_trig_idx
目的是给event_entity中的 event_trig_idx
'''
def complete_event_info(sent_event_dict, trig_dict):
    for event_idx, event_entity in sent_event_dict.items():
        main_trig = event_entity.event_trig
        cand_trig_idx = confirm_trig_idx(trig_dict, main_trig)
        assert cand_trig_idx != ""
        event_entity.event_trig_idx = cand_trig_idx

'''
思路： 过一遍 sent_event_idx,构建一个dict字典， key是trig_idx, value是event_entity
'''
def event_argu_replace(trig_dict, sent_event_dict):
    #---- 对 sent_event_dict中的每一个event.event_trig_idx赋值
    trigidx_event_dict = dict()
    for event_idx, event_entity in sent_event_dict.items():
        event_trig_idx = event_entity.event_trig_idx
        if event_trig_idx in trigidx_event_dict.keys():
            trigidx_event_dict[event_trig_idx].append(event_entity)
        else:
            trigidx_event_dict[event_trig_idx] = [event_entity]

    # --- 给论元trigger赋值 trig_idx
    for event_idx, event_entity in sent_event_dict.items():
        if isinstance(event_entity.first_argu, Trigger):  # 替换
            cand_trig_idx = confirm_trig_idx(trig_dict, event_entity.first_argu)
            if cand_trig_idx == "":
                proper_trig_entity = find_proper_trig(trig_dict, event_entity.first_argu)
                assert proper_trig_entity != None
                event_entity.first_argu = proper_trig_entity

            cand_trig_idx = confirm_trig_idx(trig_dict, event_entity.first_argu)
            event_entity.first_argu.trig_idx = cand_trig_idx

        if event_entity.second_argu != None and isinstance(event_entity.second_argu, Trigger):
            cand_trig_idx = confirm_trig_idx(trig_dict, event_entity.second_argu)
            if cand_trig_idx == "":
                proper_trig_entity = find_proper_trig(trig_dict, event_entity.second_argu)
                assert proper_trig_entity != None
                event_entity.second_argu = proper_trig_entity

            cand_trig_idx = confirm_trig_idx(trig_dict, event_entity.second_argu)
            event_entity.second_argu.trig_idx = cand_trig_idx

    # # --- 将Trigger论元不必替换成Event事件， 只需要根据arguTrig中的trig_idx, 找到对用的论元即可
    # 一些事件可能会被删除，所以新建一个dict
    for event_idx, event_entity in sent_event_dict.items():
        if isinstance(event_entity.first_argu, Protein):
            event_entity.first_argu_idx = event_entity.first_argu.prot_idx
        elif isinstance(event_entity.first_argu, Trigger): # 替换
            argutrig_idx = confirm_trig_idx(trig_dict, event_entity.first_argu)
            assert argutrig_idx != ""
            # 利用argutrig_idx 找事件
            if argutrig_idx in trigidx_event_dict.keys():
                cand_event_set = trigidx_event_dict[argutrig_idx]  # 不考虑太多，默认 cand_event_set 中有一个事件
                argu_event_idx = cand_event_set[0].event_idx
                event_entity.first_argu_idx = argu_event_idx

        if event_entity.second_argu != None and isinstance(event_entity.second_argu, Protein):
            event_entity.second_argu_idx = event_entity.second_argu.prot_idx

        elif event_entity.second_argu != None and isinstance(event_entity.second_argu, Trigger):
            argutrig_idx = confirm_trig_idx(trig_dict, event_entity.second_argu)
            assert argutrig_idx != ""
            # 利用argutrig_idx 找事件
            if argutrig_idx in trigidx_event_dict.keys():
                cand_event_set = trigidx_event_dict[argutrig_idx]  # 不考虑太多，默认 cand_event_set 中有一个事件
                argu_event_idx = cand_event_set[0].event_idx
                event_entity.second_argu_idx = argu_event_idx
    # --- 剔除掉一些不合适的event
    retain_event_dict = dict()
    for event_idx, event_entity in sent_event_dict.items():
        if event_entity.second_argu == None and event_entity.first_argu_idx != "":
            retain_event_dict[event_idx] = event_entity
        elif event_entity.second_argu != None and event_entity.second_argu_idx != "":
            retain_event_dict[event_idx] = event_entity
    return retain_event_dict
'''
argu_trig不在trig_dict 中， 要替换成trig_dict中的某个trig
argu_trig.trig_name 与trig_dict 中某个trigger的名字相同，但是位置有偏差，替换
'''
def find_proper_trig(trig_dict, argu_trig):
    argu_trig_name = argu_trig.trig_name
    for trig_idx, trig_entity in trig_dict.items():
        if trig_entity.trig_name == argu_trig_name:
            return trig_entity
    return None