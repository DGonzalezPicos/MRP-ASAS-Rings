"""
Created on Thu Apr 22 08:46:50 2021
@author: dario

Title: Panotpes client
"""
import numpy as np
import matplotlib.pyplot as plt
from panoptes_client import Panoptes, Project, SubjectSet, Subject

def panoptes_add(filelist, subject_set_name='Testing'):
    Panoptes.connect(username='DarioGonzalez', password='skyswaglord')
    
    project = Project(15697)
    # project.display_name
    # project.links.workflows
    
    subject_set = SubjectSet()
    subject_set.links.project = project
    try:
        subject_set.display_name = subject_set_name
    except:
        subject_set.display_name = subject_set_name + '-2'
    subject_set.save()
    
    new_subjects = []
    print('Merging subjects...')
    for filename in filelist:
        subject = Subject()
        
        subject.links.project = project
        subject.add_location(filename)
        subject.save()
        
        # subject.metadata.update(metadata)
        # print(filename)
    
        
        new_subjects.append(subject)
    subject_set.add(new_subjects)
    print('Saving subject set...')
    subject_set.save()
    print('Finished!!!!')
    