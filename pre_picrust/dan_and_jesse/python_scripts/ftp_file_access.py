from ftplib import FTP
from sys import argv
import os
from copy import copy
from random import shuffle

def handle_login(ftp_address):
    """Return a logged in ftp connection"""
    ftp = FTP(ftp_address)
    
    ftp.login() 
    return ftp


def change_dir(ftp_connection,new_directory):
    """Change the working directory on an ftp connection"""
    ftp_connection.cwd(new_directory)
    return ftp_connection

def download_file(ftp_connection,filename_on_server, local_filepath):
    """Download a file in the current directory from the server
    ftp_connection -- FTP object from ftplib
    
    filename_on_server -- a string for the file name (not full path) on the server.
      This needs to be in the current directory.

    local_filepath -- the local filepath to which the file should be saved
    Will overwrite without warning any existing files at this location.
    """
    f = open(local_filepath,"w+")
    
    def handle_download(block):
        f.write(block)

    ftp_connection.retrbinary('RETR %s' %(filename_on_server),callback=handle_download)
    f.close()


def find_space_needed(ftp,files_to_check,debug=False):
    """Return the final space needed to store files_to_check
    ftp -- an ftp connection, already logged in and in the correct dir
    files_to_check -- a list of filenames that are planned downloads
    """
    size_needed = 0
    
    tmp_file_list = copy(files_to_check)
    shuffle(tmp_file_list) # gives a better estimate, faster
    
    for i,file_to_check in enumerate(tmp_file_list):
         size_needed += ftp.size(file_to_check)
         if debug:
             estimate = size_needed/(i+1)*(len(tmp_file_list)-i)
             print "On file %i, estimated total: %s" %(i+1,bytes_to_size_string(size_needed + estimate))
    return size_needed     


def bytes_to_size_string(bytes_needed,kb_size=1024,\
        unit_abbreviations ={"byte":"b","kilobyte":"Kb",\
          "megabyte":"MB","gigabyte":"GB",\
          "terabyte":"TB"}):
    """Return a descriptive string from a number of bytes
    kb_size -- integer for the number of bytes in a kb (either 1024 or 1000)
    unit_abbreivations -- dict mapping unit names (e.g. 'kilobyte') to 
    the desired abbreviation (e.g. 'Kb)
    """
    
    kilobyte = kb_size 
    megabyte = kilobyte**2
    gigabyte = kilobyte**3
    terabyte = kilobyte**4
     
    if bytes_needed < kilobyte:
            unit = "byte"
            scalar = bytes_needed
    elif bytes_needed < megabyte:
            unit = "kilobyte"
            scalar = bytes_needed/kilobyte
    elif bytes_needed < gigabyte:
            unit = "megabyte"
            scalar = bytes_needed/megabyte
    elif bytes_needed < terabyte:
            unit = "gigabyte"
            scalar = bytes_needed/gigabyte
    else:
            unit = "terabyte" 
            scalar = bytes_needed/terabyte

    return  "%i %s" %(scalar,unit_abbreviations.get(unit))


def login_and_download_files(ftp_site,ftp_folder,file_names,output_dir='./',debug=False):
    """Login anonymously to an ftp site and download a list of specified files
    ftp_site -- ftp address of the site (e.g. 'ftp.jgi-psf.org')
    ftp_folder -- ftp folder of the files (e.g. '/pub/IMG/img_w_v340')
    file_names -- a *list* of files to download (e.g. ['00.taxon.tab.txt'])
    local_filepath -- dir to save files (e.g. './')
    """
    
    if not os.path.isdir(output_dir):
        raise RuntimeError("Specified output_dir '%s' does not appear to be a directory." % output_dir)
    
    if debug:
        print "Logging in to: %s" %ftp_site
    ftp = handle_login(ftp_site)
    
    if debug:
        print "Changing to the '%s' directory..." %ftp_folder
    
    ftp = change_dir(ftp, ftp_folder)
    
    if debug:
        print "Checking file sizes..."
        bytes_needed = find_space_needed(ftp,file_names)
        size_description = bytes_to_size_string(bytes_needed)
        print "%s needed to download all %i files" %(size_description,len(file_names))
        

    for i,file_to_download in enumerate(file_names):
        local_filepath = os.path.join(output_dir, file_to_download)
        if debug:
            print "Downloading %s (file %i/%i) from %s to: %s" %(file_to_download,\
              i+1,len(file_names),ftp_site, local_filepath)
        
        download_file(ftp, file_to_download, local_filepath)      

    #Politely close connection
    ftp.quit() 



def parse_IMG_index_file(lines):
    """Parse img lines yielding a fields"""
    header_line_index = 0
    for i,line in enumerate(lines):
        if i == header_line_index:
            continue
        fields = line.strip().split("\t")
        taxon_oid,taxon_display_name,domain,seq_status = fields
        yield taxon_oid



def taxon_ids_from_ftp_index_file(ftp_site = 'ftp.jgi-psf.org',\
                                 ftp_folder = '/pub/IMG/img_w_v340',\
                                 index_file='00.taxon.tab.txt',
                                 index_file_parser=parse_IMG_index_file,\
                                 output_dir="./",debug=False):
    """Download and parse a list of taxon ids from an index file.
    tmp_dir = the directory to store IMG's readme file)
    """
    
    # Download index file
    login_and_download_files(ftp_site=ftp_site,ftp_folder=ftp_folder,\
      file_names=[index_file],\
      output_dir = output_dir)
    
    # Parse index file to retrieve taxon ids
    f = open(os.path.join(output_dir,index_file),"U")
    taxon_id_generator = index_file_parser(f.readlines())
    taxon_oids = [id for id in taxon_id_generator]
    
    return taxon_oids


def download_img_genomes(ftp_site = 'ftp.jgi-psf.org',\
    ftp_folder = '/pub/IMG/img_w_v340',\
    index_file = '00.taxon.tab.txt',\
    index_file_parser = parse_IMG_index_file,\
    output_dir = './',limit = None, debug = False):

    # Download and parse the list of available genomes
    if debug:
        print "Downloading a list of available genomes...."

    taxon_ids_to_download =\
      taxon_ids_from_ftp_index_file(ftp_site=ftp_site,\
        ftp_folder = ftp_folder, index_file = index_file,\
        output_dir = output_dir, debug = debug)

    if debug:
        print "Found %i available genomes." % len(taxon_ids_to_download)
        print "Starting download...."
    
    if limit:
        taxon_ids_to_download = \
          taxon_ids_to_download[0:max(limit,1)]
        if debug:
            print "Download limited to %i genomes" % limit

    download_img_genomes_by_taxonid(taxon_ids_to_download,\
      output_dir = output_dir, debug=debug)

def download_img_genomes_by_taxonid(taxon_ids,output_dir='./',debug = False):
    """Download an img genome by taxon_id
    taxon_id -- the taxon_id of the genome to download
    outpath -- the path to which the genome .tar.gz file should be saved
               if set to None, defaults to the taxonid + .tar.gz in the current folder
    """
    files_to_download = ["%s.tar.gz" % str(taxon_id) for taxon_id in taxon_ids]
    
    login_and_download_files(ftp_site='ftp.jgi-psf.org',\
      ftp_folder = '/pub/IMG/img_w_v340',\
      file_names = files_to_download,\
      output_dir = './',\
      debug = debug)
     

if __name__ == "__main__":
    # For prototype purposes, let's grab an IMG genome
    print "Example of downloading an individual genome by taxon id"
    m_genitalium_g37_taxon_id = 638341133
    download_img_genomes_by_taxonid([m_genitalium_g37_taxon_id],"./m_genitalium.tar.gz",debug=True)
    
    # Get the IMG genome list
    print "Retrieving taxon oids..."
    download_img_genomes(debug=True) 
