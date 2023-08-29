from entity import Protein, Event
def data2json(prot_dict, trig_dict, event_dict, sent_json_form):
    entity_mention_set = []
    for prot_idx, prot_entity in prot_dict.items():
        prot_mentions = {}
        prot_mentions["id"] = prot_idx
        prot_mentions["text"] = prot_entity.prot_name
        prot_mentions["entity_type"] = "Protein"
        prot_mentions["start"] = prot_entity.prot_start
        prot_mentions["end"] = prot_entity.prot_end+1
        assert len(prot_mentions["text"].split()) == prot_mentions["end"]-prot_mentions["start"]
        entity_mention_set.append(prot_mentions)

    for trig_idx, trig_entity in trig_dict.items():
        trig_mentions = {}
        trig_mentions["id"] = trig_idx
        trig_mentions["text"] =trig_entity.trig_name
        trig_mentions["entity_type"] =trig_entity.trig_type
        trig_mentions["start"] =trig_entity.trig_start
        trig_mentions["end"] = trig_entity.trig_end+1
        assert len(trig_mentions["text"].split()) == trig_mentions["end"]-trig_mentions["start"]
        entity_mention_set.append(trig_mentions)
    sent_json_form['entity_mentions'] = entity_mention_set

    event_mention_set = []
    for event_idx, event_entity in event_dict.items():
        event_mentions = {}
        event_mentions["id"] = event_idx
        event_mentions["event_type"] = event_entity.event_type
        event_trig = event_entity.event_trig
        event_mentions["trigger"] = {"text": event_trig.trig_name, "start": event_trig.trig_start, "end": event_trig.trig_end+1}

        argument_set = []
        first_argu = event_entity.first_argu
        first_argu_type = event_entity.first_argu_type
        get_argument_info(first_argu, first_argu_type, argument_set)

        second_argu = event_entity.second_argu
        second_argu_type = event_entity.second_argu_type
        if second_argu != None:
            get_argument_info(second_argu, second_argu_type, argument_set)

        event_mentions["arguments"] = argument_set
        event_mention_set.append(event_mentions)
    sent_json_form["event_mentions"] = event_mention_set
#
def get_argument_info(argument, argu_type, argument_set):
    argument_info = {}
    if isinstance(argument, Protein):
        argument_info["entity_id"] = argument.prot_idx
        argument_info["text"] = argument.prot_name
        argument_info["role"] = argu_type
        argument_set.append(argument_info)
    elif isinstance(argument, Event):
        event_trig = argument.event_trig
        argument_info["entity_id"] = event_trig.trig_idx
        argument_info["text"] = event_trig.trig_name
        argument_info["role"] = argu_type
        argument_set.append(argument_info)
    else:
        print("what else, json_utils.py第66行")

# test 数据集
def testdata2json(prot_dict, sent_json_form):
    entity_mention_set = []
    for prot_idx, prot_entity in prot_dict.items():
        prot_mentions = {}
        prot_mentions["id"] = prot_idx
        prot_mentions["text"] = prot_entity.prot_name
        prot_mentions["entity_type"] = "Protein"
        prot_mentions["start"] = prot_entity.prot_start
        prot_mentions["end"] = prot_entity.prot_end + 1
        assert len(prot_mentions["text"].split()) == prot_mentions["end"] - prot_mentions["start"]
        entity_mention_set.append(prot_mentions)

    sent_json_form['entity_mentions'] = entity_mention_set
    sent_json_form["event_mentions"] = []


