"""Utilities to download data from NCBI SRA"""
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)
import sys
from subprocess import Popen, PIPE, STDOUT


def run_download_sra_script(download_root_location=None,
                            ascp_key_path=None,
                            srp_id_file=None,
                            srp_id_list=None):
    """Download data from SRA.

    Parameters
    ------------
    download_root_location : string
                             Path to download SRA files
    ascp_key_path : string
                    Location for aspera private keypp
    srp_id_list : list
                  List of SRP ids for download

    srp_id_file : string
                  File containing list of SRP Ids, one per line

    """
    cmd = 'download_sra_data'
    if download_root_location:
        cmd += ' --out {} '.format(download_root_location)
    if ascp_key_path:
        cmd += ' --ascp {}'.format(ascp_key_path)
    if srp_id_file:
        cmd += ' --file {}'.format(srp_id_file)
    elif srp_id_list:
        cmd += ' '
        cmd += ' '.join(srp_id_list)
    cmds = cmd.strip().split(' ')
    proc = Popen(cmds, stdout=PIPE, stderr=STDOUT)
    while True:
        output = proc.stdout.readline()
        if output == '' and proc.poll() is not None:
            break
        if output:
            print(str(output.strip(), 'utf-8'))
    rc = proc.poll()
    return rc
