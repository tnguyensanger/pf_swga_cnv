'''
For each of the PF SWGA X10 samples, reimport the WGS samples.  The old WGS samples were lost because lustre scratch109 was decommissioned.
'''
import os
import subprocess
import csv
import json
from collections import namedtuple
import argparse
import logging
import sys



TARGET_CRAM_FOFN = os.path.dirname(os.path.realpath(__file__)) + os.sep + os.pardir + os.sep + os.pardir + os.sep + "data" + os.sep + "pf_swga_cnv_wgs_import.in.fofn"
SAMPLE_MAP_CSV = os.path.dirname(os.path.realpath(__file__)) + os.sep + os.pardir + os.sep + os.pardir + os.sep + "data" + os.sep + "sWGA_plexing_sample_map.csv"

LOGGER = logging.getLogger(__name__)
ch = logging.StreamHandler(sys.stdout)
ch.setLevel(logging.DEBUG)
# create formatter and add it to the handlers
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
ch.setFormatter(formatter)
# add the handlers to the logger
LOGGER.addHandler(ch)
LOGGER.setLevel(logging.DEBUG)
# root = logging.getLogger()
# root.addHandler(ch)
# root.setLevel(logging.DEBUG)


RunInfo = namedtuple("RunInfo", ["run", "lane_index", "tag_index"])

class IrodsBaton:
    """
    Uses Irods Baton program which uses jq which pipes json to baton-metaquery
    """
    # jq -n '{collection: "/humgen/projects/crohns/20130909", avus: [{attribute: "study_internal_id", value: "2323"}, {attribute: "target", value: "1", o: "="}]}' | /software/solexa/pkg/baton/current/bin/baton-metaquery --avu

    
    @staticmethod
    def list_sample(sample=None):
        """
        Lists the samples in alphabetical order that were included in the run.

           
        """
        
        cmd = ["baton-metaquery", "--avu", "--obj"]
        baton_input = {"collection": "/seq", 
                       "avus": [{"attribute": "sample", "value": sample}, 
                                {"attribute": "manual_qc", "value": "1", "o": "="}, 
                                {"attribute": "target", "value": "1", "o": "="},
                                {"attribute": "type", "value": "cram"}
                                ] 
                       }

            
        baton_input_str = json.dumps(baton_input)
        LOGGER.debug(baton_input_str)
        jqproc = subprocess.Popen(["jq", "-n", baton_input_str], stdout=subprocess.PIPE)
        baton_output = subprocess.check_output(cmd, stdin=jqproc.stdout)

        objects_list = json.loads(baton_output)
        if len(objects_list) > 1:
            LOGGER.warn(objects_list)
            raise ValueError("There are more than 1 files for that sample")
        
        return objects_list[0]
    
    



def create_import_fofn(cram_output_dir, output_fofn):
    with open(SAMPLE_MAP_CSV, 'r') as fh_in, open(output_fofn, 'w') as fh_out:
        reader = csv.DictReader(fh_in)
        writer = csv.DictWriter(fh_out, delimiter="\t", fieldnames=["path", "sample", "library", "lane", "study"])
        writer.writeheader()
        for inrow in reader:
            wgs_sample_name = inrow["seqscape_name_wgs_Hiseq2000"]

        
                        # EG of object_list:
#                 [{u'avus': [{u'attribute': u'type', u'value': u'cram'}, 
#                             {u'attribute': u'sample_id', u'value': u'2408537'}, 
#                             {u'attribute': u'library', u'value': u'15108295'}, 
#                             {u'attribute': u'md5', u'value': u'b54083db43f7b290de69daf26b2a53ce'}, 
#                             {u'attribute': u'is_paired_read', u'value': u'1'}, 
#                             {u'attribute': u'ebi_sub_date_history', u'value': u'[2016-01-22T01:11:28] 2016-01-19'},
#                              {u'attribute': u'study', u'value': u'IHTP_PWGS 1151-PF-GH-AMENGA-ETEGO'}, 
#                              {u'attribute': u'total_reads', u'value': u'22254088'}, 
#                              {u'attribute': u'ebi_run_acc', u'value': u'ERR1214170'}, 
#                              {u'attribute': u'manual_qc', u'value': u'1'}, 
#                              {u'attribute': u'sample_donor_id', u'value': u'3909STDY6204252'}, 
#                              {u'attribute': u'library_type', u'value': u'Standard'}, 
#                              {u'attribute': u'sample_common_name', u'value': u'Plasmodium falciparum'}, 
#                              {u'attribute': u'id_run', u'value': u'18091'}, 
#                              {u'attribute': u'ebi_sub_acc', u'value': u'ERA553489'}, 
#                              {u'attribute': u'tag_index', u'value': u'8'},
#                               {u'attribute': u'study_accession_number', u'value': u'ERP000190'}, 
#                               {u'attribute': u'reference', u'value': u'/lustre/scratch110/srpipe/references/Plasmodium_falciparum/3D7/all/bwa/MAL.version2.1.4.fasta'}, 
#                               {u'attribute': u'study_id', u'value': u'3909'}, 
#                               {u'attribute': u'lane', u'value': u'4'}, 
#                               {u'attribute': u'target', u'value': u'1'}, 
#                               {u'attribute': u'study_title', u'value': u'MALARIA SURVEILLANCE, GENOMIC EPIDEMIOLOGY AND INVESTIGATION OF HOST-PARASITE INTERACTION'}, 
#                               {u'attribute': u'sample', u'value': u'3909STDY6204252'}, 
#                               {u'attribute': u'sample_accession_number', u'value': u'ERS902136'}, 
#                               {u'attribute': u'ebi_sub_md5', u'value': u'b54083db43f7b290de69daf26b2a53ce'}, 
#                               {u'attribute': u'alignment', u'value': u'1'}, 
#                               {u'attribute': u'library_id', u'value': u'15108295'}, 
#                               {u'attribute': u'ebi_sub_date', u'value': u'2016-01-21'}, 
#                               {u'attribute': u'sample_supplier_name', u'value': u'PF1195-Cx'}], 
#                   u'data_object': u'18091_4#8.cram', 
#                   u'collection': u'/seq/18091'}]
            object_list = IrodsBaton.list_sample(sample=wgs_sample_name)
            LOGGER.debug(object_list)
            outrow = {}
            for obj_detail, obj_det_val in object_list.items():
                if obj_detail == "data_object": 
                    path_basename = obj_det_val
                if obj_detail == "avus":
                    for attr_dict in obj_det_val:
                        attribute = attr_dict["attribute"]
                        value = attr_dict["value"]
                        if attribute == "library":
                            outrow["library"] = value
                        elif attribute == "study":
                            outrow["study"] = value
                        elif attribute == "sample":
                            outrow["sample"] = value

                    
            outrow["path"] = cram_output_dir  + os.sep + path_basename
            outrow["lane"] = path_basename.replace(".cram", "")
            writer.writerow(outrow)

                    

                
                
if __name__ == '__main__':
    
    # http://stackoverflow.com/questions/3991104/very-large-input-and-piping-using-subprocess-popen

    parser = argparse.ArgumentParser()
    parser.add_argument("-cram_output_dir", help="output directory for the crams")
    parser.add_argument("-output_fofn", help="output FOFN with metadata to write to")
    parser.description = "Create a fofn with metadata you can use for import the corresponding WGS Hiseq2000 samples."


    args = parser.parse_args()

    create_import_fofn(cram_output_dir=args.cram_output_dir, output_fofn=args.output_fofn)
   
    
    LOGGER.info("Done!")
                
           
                
        
            