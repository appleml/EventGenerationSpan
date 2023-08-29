# !/usr/bin/env python
# -*- coding:utf-8 -*-
'''
处理数据： 与Text2Event同样的树状结构数据， 应对重复的实体，就用临近原则
处理步骤： (1) data_json.py (2) convert_text_to_tree.py
'''
import os
import json
from collections import Counter, defaultdict
from data_convert.format.text2tree import Text2Tree
from data_convert.task_format.event_extraction import Event
from data_convert.utils import read_file, check_output, data_counter_to_table, get_schema, output_schema
from nltk.corpus import stopwords

english_stopwords = set(stopwords.words('english') + ["'s", "'re", "%"])

def convert_file_tuple(file_tuple, data_class=Event, target_class=Text2Tree,
                       output_folder='data/2011_tree_data/',
                       ignore_nonevent=False, zh=False,
                       mark_tree=False, type_format='subtype'):
    counter = defaultdict(Counter)
    data_counter = defaultdict(Counter)

    event_schema_set = set()

    span_output_folder = output_folder + '_span'
    if not os.path.exists(span_output_folder):
        os.makedirs(span_output_folder)

    for in_filename, output_filename in file_tuple(output_folder):
        span_output_filename = output_filename.replace(output_folder, span_output_folder)
        event_output = open(output_filename + '.json', 'w')
        span_event_output = open(span_output_filename + '.json', 'w')

        for line in read_file(in_filename):
            document = data_class(json.loads(line.strip()))
            # if document.doc_key == "PMID-8496329.txt":
            #     print("############################")
            for sentence in document.generate_sentence(type_format=type_format):

                if ignore_nonevent and len(sentence['events']) == 0:
                    continue
                # source是句子， target是生成的事件字符串
                source, target = target_class.annotate_predicate_arguments(
                    tokens=sentence['tokens'],
                    predicate_arguments=sentence['events'],
                    zh=zh
                )

                for event in sentence['events']:
                    event_schema_set = event_schema_set | get_schema(event)
                    sep = '' if zh else ' '
                    predicate = sep.join([sentence['tokens'][index] for index in event['tokens']])
                    counter['pred'].update([predicate])
                    counter['type'].update([event['type']])
                    data_counter[in_filename].update(['event'])
                    for argument in event['arguments']:
                        data_counter[in_filename].update(['argument'])
                        counter['role'].update([argument[0]])

                data_counter[in_filename].update(['sentence'])

                event_output.write(json.dumps({'text': source, 'event': target}, ensure_ascii=False) + '\n')

                span_source, span_target = target_class.annotate_span(
                    tokens=sentence['tokens'],
                    predicate_arguments=sentence['events'],
                    zh=zh,
                    mark_tree=mark_tree
                )

                span_event_output.write(json.dumps({'text': span_source, 'event': span_target}, ensure_ascii=False) + '\n')

        event_output.close()
        span_event_output.close()

        check_output(output_filename)
        check_output(span_output_filename)
        print('\n')
    output_schema(event_schema_set, output_file=os.path.join(output_folder, 'event.schema'))
    output_schema(event_schema_set, output_file=os.path.join(span_output_folder, 'event.schema'))
    print('Pred:', len(counter['pred']), counter['pred'].most_common(10))
    print('Type:', len(counter['type']), counter['type'].most_common(10))
    print('Role:', len(counter['role']), counter['role'].most_common(10))
    print(data_counter_to_table(data_counter))
    print('\n\n\n')


def convert_bio_event(output_folder='data/2011_tree_data/', type_format='subtype', ignore_nonevent=False, mark_tree=False):
    from data_convert.task_format.event_extraction import bioevent_file_tuple
    convert_file_tuple(file_tuple=bioevent_file_tuple,
                       output_folder=output_folder,
                       ignore_nonevent=ignore_nonevent,
                       mark_tree=mark_tree,
                       type_format=type_format,
                       )

if __name__ == "__main__":
    absolute_path = "/home/fangsu/myworks/EventGeneration/"
    #
    # train_data = absolute_path + "data/2011_data_one/2011_train"
    # train_genia = absolute_path + "data/2011_genia_parse_one_split_bio/2011_train"
    # train_processed = absolute_path + "data/2011_tree_data/2011_train"  # 转换为增广后的数据
    # read_data(train_data, train_genia, train_processed)
    #
    # devel_data = absolute_path + "data/2011_data_one/2011_devel"
    # devel_genia = absolute_path + "data/2011_genia_parse_one_split_bio/2011_devel"
    # devel_processed = absolute_path + "data/2011_tree_data/2011_devel"
    # read_data(devel_data, devel_genia, devel_processed)
    #
    # test_data = absolute_path + "data/2011_data_one/2011_test"
    # test_genia = absolute_path + "data/2011_genia_parse_one_split_bio/2011_test"
    # test_processed = absolute_path + "data/2011_tree_data/2011_test"
    # read_test_data(test_data, test_genia, test_processed)


    type_format_name = 'subtype'
    convert_bio_event(absolute_path + "data/2011_tree_%s" % type_format_name,
                          type_format=type_format_name,
                          ignore_nonevent=False,
                          mark_tree=False
                          )

    convert_bio_event(absolute_path + "data/2011_tree_%s" % type_format_name,
                          type_format=type_format_name,
                          ignore_nonevent=False,
                          mark_tree=False
                          )

    convert_bio_event(absolute_path + "data/2011_tree_%s" % type_format_name,
                          type_format=type_format_name,
                          ignore_nonevent=False,
                          mark_tree=False
                          )