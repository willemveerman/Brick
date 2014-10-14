# -*- coding: utf-8 -*-
"""
@author: willem
"""

import urllib2
from Bio.Seq import Seq
from bioservices import *
import xml.etree.ElementTree as ET

class Brick():
    """Interface to the iGEM Parts Repository - http://parts.igem.org/Registry_API

        Enter part ID in quotes, e.g.:
        
        Brick("BBa_C0079")

    """
    def __init__(self,id):  
        self.id = id
        try:
            #Access part URL and extract XML
            biobrick_xml_raw = urllib2.urlopen("http://parts.igem.org/cgi/xml/part.cgi?part="+self.id)
            
            #Parse XML and get root with ElementTree
            self.parsed = ET.parse(biobrick_xml_raw)
            self.root = self.parsed.getroot()
            
            #Check that part exists, if so print message. If not, pass error
            for part in self.root.findall('part_list'):
                if part.find('ERROR') is not None:
                    print part.find('ERROR').text
                else:
                    print "Part successfully loaded."                
        except urllib2.URLError:
            print "No connection to the internet at present."
               
    def __str__(self):
        return "BioBrick object, ID: {} - for further information, use 'overview' method.".format(self.id)
        
    def part_attrib(self,node):
        """ This function retrieves the part's attributes from its XML file.
        
            Attributes which are directly under the 'part_list' tree can
            be called by simply typing their name in quotes, e.g.:
            
            Brick.part_attrib('part_type')
            
            The attributes which can be accessed in this way are:
            
            part_id
            part_name
            part_short_name
            part_short_desc
            part_type
            uniprot_id
            release_status
            sample_status
            part_results
            part_nickname
            part_rating
            part_url
            part_entered
            part_author
            deep_subparts
            specified_subparts
            specified_subscars
            
            Other parts can be accessed via their 'address', i.e. 'sequences/seq_data'
            
        """
        
        if node == 'uniprot_id':
            swisspro_value = None
            for parameter in self.root.iter(tag='parameter'):
                name = parameter.find('name')
                if name is not None and name.text == 'swisspro':
                    swisspro_value = parameter.find('value').text
                    break
                
            return swisspro_value or "No UniProt ID present."
        #The block above allows the user to enter the query 'uniprot_id',
        #despite the fact that that node is not on the correct tree
        #If uniprot_id is entered, the code finds the node 'swisspro' and enters
        #its corresponding value, found on the node directly underneath
        #If there is no 'swisspro' node then an error message is returned.
        
        else:
            if [part.find('{}'.format(node)) for part in self.root.find('part_list')][0] is not None:
                return [part.find('{}'.format(node)).text for part in self.root.find('part_list')][0].replace("\n","")
            else:
                return "No {} node present in XML.".format(node)
        #The code above searches for the node enterted as the function's argument and returns the node's value
        #If the node is not found, an error message is returned
        
    def seq(self):
        """For convenience, the method call which outputs the part's sequence;
        
            Brick.part_attrib(self,'sequences/seq_data')
            
            is here bundled as a method in itself.
        """
        return Brick.part_attrib(self,'sequences/seq_data')
        
    def overview(self):
        """Quickly and conveniently prints the most relevant information
            
            by calling Brick.part_attrib() method several times.
        
        """
        print "Part ID  :",Brick.part_attrib(self,'part_name')
        print "Part Type:",Brick.part_attrib(self,'part_type')
        print "Part Nick:",Brick.part_attrib(self,'part_nickname')
        print "Part Desc:",Brick.part_attrib(self,'part_short_desc')
        print "Part URL :",Brick.part_attrib(self,'part_url')
        print "Part Seq :\n",Brick.part_attrib(self,'sequences/seq_data')


class Protein(Brick):
    """This class is a sub-class of Brick() specifically for coding sequences.

        As with Brick(), enter part ID in quotes, e.g.:
        
        Protein("BBa_C0079")
        
        Non-coding sequences can be entered as Protein() objects
        
        However, an error will be returned if they have no UniProt ID in their XML.

    """
    def __init__(self,id):
        Brick.__init__(self,id)
     
    def __str__(self):
        return "BioBrick protein object, ID: {} - for further information, use 'overview' method.".format(self.id)
        
    def protein_seq(self):
        """Part sequence translated to protein via BioPython.
        
            Returns a SeqRecord object.
        """
        return Seq(Brick.seq(self)).translate()
        
    def uniprot_overview(self,format):
        """Uses BioServices to query UniProt for a UniProt Id.
        
            format paramter - 'txt', 'xml', 'rdf', 'gff' or 'fasta'
        """
        if Brick.part_attrib(self,'uniprot_id') == "No UniProt ID present.":
            return "No UniProt ID present."
        else:
            #Initialize BioSerives UniProt object
            bioservices_up_obj = UniProt()
            
            #Search UniProt by part's UniProt ID
            results = bioservices_up_obj.searchUniProtId(str(Brick.part_attrib(self,'uniprot_id')),format)
            return results
            
    def structures(self,c=None,detail=None,format=None):
        """Uses BioServices to find PDB accession IDs and files for a UniProt ID.
            
            For a list of structure names, omit all arguments.            
            
            For number of structures, set param c='count'.
            
            For detail of a particular structure: 
            1. specify structure number with param detail (1-indexed)
            2. specify format - 'FASTA', 'pdb', 'cif' or 'xml'
            
            For example, 1st structure in XML:
            
            Protein.structures(detail=1,format='xml')
            
        """
        if Brick.part_attrib(self,'uniprot_id') == "No UniProt ID present.":
            return "No structures yet ascertained."
        else:
            #Initialize BioServices UniProt object
            bioservices_up_obj = UniProt()
            
            #Use mapping method to find PDB IDs for this part's UniProt ID
            #Ouput is a dictionary with UniProt ID as 1st key and list of PDB IDs as its value
            results = bioservices_up_obj.mapping("ID","PDB_ID",str(Brick.part_attrib(self,'uniprot_id')))
            
            if c is 'count':
                #Output number of structures - len of list which is 1st value of results dict
                return len(results[str(Brick.part_attrib(self,'uniprot_id'))])
            elif detail is not None and format is not None:
                
                #Initialize BioServices PDB object
                bioservices_pdb_obj = PDB()
                #Use getFile method to retrieve PDB file - chooses one list element
                return bioservices_pdb_obj.getFile(results[str(Brick.part_attrib(self,'uniprot_id'))][detail-1],format)
            else:
                #Return dict values; PDB IDs
                return results.values()
        
    def go_attributes(self):
        """Uses BioServices to query QuickGO Gene Ontology database.
        
            Returns a Pandas.DataFrame, displaying only columns 4 & 5.
            
        """
        #Initialize BioServices QuickGO object
        bioservices_quickgo_obj = QuickGO()
        #Search QuickGO for protein UniProt ID
        res = bioservices_quickgo_obj.Annotation_from_protein(protein=str(Brick.part_attrib(self,'uniprot_id')))
        
        #Use Pandas.DataFrame method object iloc to select specific columns
        print res.iloc[:,[4,5]]
        
    def get_models(self,x=None,lim=None):
        """Uses BioServices to search BioModels for models that contain proteins
            which share a specified GO-ID with the part contained in the object.
        
            First, call the Protein.go_attributes() method.
            Then, choose a row GO-ID to query.
            Specify this row with the x parameter.
            Specify the maximum number of proteins from UniProt containing the GO-ID
            with lim.
            
            For example:
            
            Protein.get_models(x=1,lim=100)
            
            Will return all models which contain any of the 1st 100 proteins in UniProt
            which are associated with GO-ID 1 from Protein.go_attributes().
            
        """
        #Initialize BioServices objects
        bioservices_up_obj = UniProt()
        bioservices_quickgo_obj = QuickGO()
        bioservices_biomodels_obj = BioModels()
        
        #Search QuickGO for UniProt ID of part
        res = bioservices_quickgo_obj.Annotation_from_protein(protein=str(Brick.part_attrib(self,'uniprot_id')))
        
        #Prepare container for part GO-IDs
        go_id = []
        
        #Output number of GO-IDs
        go_number = len(res['goID'])
        
        #Append all GO-IDs to a list
        for i in range(go_number):
            go_id.append(str(res.iloc[i]['goID']))
        
        #Search UniProt for one item from list. Output is dict.
        results = bioservices_up_obj.quick_search(query=str(go_id[x]),limit=lim)
        
        #Remove the part itself from list of UniProt results.
        for i in results.keys():
            if i == Brick.part_attrib(self,'uniprot_id'):
                del results[i]
                
        #Create variables to store results
        model_list = []
        protein_list =[]
        
        #For each protein that is associated with the GO-ID,
        #search for it in BioModels.
        #If it's found, append the model ID to one list and protein ID to another.
        for i in results.keys():
            if bioservices_biomodels_obj.getModelsIdByUniprot(i):
                model_list.append(bioservices_biomodels_obj.getModelsIdByUniprot(i))
                protein_list.append(i)
           
        return model_list, protein_list