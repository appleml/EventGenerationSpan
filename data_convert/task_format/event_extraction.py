#!/usr/bin/env python
# -*- coding:utf-8 -*-
import os
from data_convert.task_format.task_format import TaskFormat

class Event(TaskFormat):
    """
    {
        "doc_id": "NYT_ENG_20130914.0094",
        "sent_id": "NYT_ENG_20130914.0094-1",
        "tokens": ["LARGO", "\u2014", "A", "judge", "on", "Friday", "refused", "to", "stop", "''", "Hiccup", "Girl", "\u2019'", "Jennifer", "Mee", "from", "giving", "interviews", "to", "the", "media", "in", "the", "final", "days", "before", "her", "murder", "trial", "."],
        "entities": [
            {"entity_id": "NYT_ENG_20130914.0094-1-8-42", "entity_type": "PER", "mention_type": "NOM", "start": 3, "end": 4, "text": "judge"},
            {"entity_id": "NYT_ENG_20130914.0094-1-39-2096", "entity_type": "PER", "mention_type": "NAM", "start": 10, "end": 12, "text": "Hiccup Girl"},
            {"entity_id": "NYT_ENG_20130914.0094-1-39-48", "entity_type": "PER", "mention_type": "NAM", "start": 13, "end": 15, "text": "Jennifer Mee"},
            {"entity_id": "NYT_ENG_20130914.0094-1-39-54", "entity_type": "PER", "mention_type": "PRO", "start": 26, "end": 27, "text": "her"}
        ],
        "relations": [],
        "events": [
            {
                "event_id": "NYT_ENG_20130914.0094-1-25-998", "event_type": "contact", "event_subtype": "broadcast",
                "trigger": {"text": "interviews", "start": 17, "end": 18},
                "arguments": [
                    {"entity_id": "NYT_ENG_20130914.0094-1-39-48", "role": "entity", "text": "Jennifer Mee"}
                ]
            },
            {
                "event_id": "NYT_ENG_20130914.0094-1-1040-1019", "event_type": "justice", "event_subtype": "trialhearing",
                "trigger": {"text": "trial", "start": 28, "end": 29},
                 "arguments": [{"entity_id": "NYT_ENG_20130914.0094-1-39-54", "role": "defendant", "text": "her"}
                 ]
            }
        ],
        "start": 231, "end": 380,
         "text": "LARGO \u2014 A judge on Friday refused to stop ''Hiccup Girl\u2019' Jennifer Mee from giving interviews to the media in the final days before her murder trial."
    }
    """

    def __init__(self, doc_json): # doc_json是dict()
        self.doc_key = doc_json['doc_id']
        self.sentence = doc_json['tokens']
        self.entities = {entity['id']: entity for entity in doc_json['entity_mentions']}
        #self.relations = doc_json['relation_mentions']
        self.events = doc_json['event_mentions']

    def generate_sentence(self, type_format='subtype'):
        events = list()
        for event in self.events:
            arguments = list()
            for argument in event['arguments']:
                argument_entity = self.entities[argument['entity_id']]
                arguments += [[argument['role'], list(range(argument_entity['start'], argument_entity['end']))]]

            event_type = event['event_type']
            events += [{
                'type': event_type,
                'tokens': list(range(event['trigger']['start'], event['trigger']['end'])),
                'arguments': arguments
            }]

        yield {'tokens': self.sentence, 'events': events}


def bioevent_file_tuple(output_folder):
    if not os.path.exists(output_folder):
        os.makedirs(output_folder, exist_ok=True)

    input_folder = "/home/fangsu/myworks/EventGeneration/data/2011_json_data"

    file_tuple = [
        (input_folder + "/train.json", output_folder + '/train'),
        (input_folder + "/devel.json", output_folder + '/devel'),
        (input_folder + "/test.json", output_folder + '/test'),
    ]

    return file_tuple


if __name__ == "__main__":
    pass
