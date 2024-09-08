# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" A collection of functions for running hmmsearch.
"""

import logging
import os
from typing import List

from helperlibs.wrappers.io import TemporaryDirectory

from .base import execute, get_config, SearchIO

hmmsearch_time_1=0
pyhmmsearch_time_1=0
hmmsearch_time_2=0
pyhmmsearch_time_2=0

def run_hmmsearch(query_hmmfile: str, target_sequence: str, use_tempfile: bool = False
                  ) -> List[SearchIO._model.query.QueryResult]:  # pylint: disable=protected-access
    """ Run hmmsearch on a HMM file and a fasta input

        Arguments:
            query_hmmfile: the path to the HMM file
            target_sequence: the fasta input to search as a string
            use_tempfile: if True, a tempfile will be written for the fasta input
                          instead of piping

        Returns:
            a list of hmmsearch results as parsed by SearchIO
    """
    import time
    from pyhmmer import easel, hmmer, plan7

    global hmmsearch_time_1
    global pyhmmsearch_time_1
    global hmmsearch_time_2
    global pyhmmsearch_time_2


    config = get_config()
    command = [config.executables.hmmsearch, "--cpu", str(config.cpus),
               "-o", os.devnull,  # throw away the verbose output
               "--domtblout", "result.domtab",
               query_hmmfile]

    # Allow for disabling multithreading for HMMer3 calls in the command line
    if config.get('hmmer3') and 'multithreading' in config.hmmer3 and \
            not config.hmmer3.multithreading:
        command = command[0:1] + command[3:]

    with TemporaryDirectory(change=True):
        try:
            if use_tempfile:
                with open("input.fa", "w", encoding="utf-8") as handle:
                    handle.write(target_sequence)
                command.append("input.fa")

                start_time = time.time()
                run_result = execute(command)
                search_time = (time.time() - start_time)
                hmmsearch_time_1 += search_time
                logging.info(f"hmmsearch_time_1: {hmmsearch_time_1}")


                # pyhmmer search begins
                start_time = time.time()
            
                hmm_file = plan7.HMMFile(query_hmmfile)
                seq_file = easel.SequenceFile("input.fa", digital=True, alphabet=easel.Alphabet.amino())
                hits_list = list(hmmer.hmmsearch(hmm_file, seq_file, cpus=80))

                pyhmmer_search_fp = open("result_pyhmmer.domtab", "wb")
                hits_list[0].write(pyhmmer_search_fp, format="domains", header=True)
                for hit in hits_list[1:]:
                    hit.write(pyhmmer_search_fp, format="domains", header=False)
                pyhmmer_search_fp.close()

                search_time = (time.time() - start_time)
                pyhmmsearch_time_1 += search_time
                logging.info(f"pyhmmsearch_time_1: {pyhmmsearch_time_1}")
                # pyhmmer search ends
            else:
                command.append('-')

                start_time = time.time()                
                run_result = execute(command, stdin=target_sequence)
                search_time = (time.time() - start_time)
                hmmsearch_time_2 += search_time
                logging.info(f"hmmsearch_time_2: {hmmsearch_time_2}")

                # pyhmmer search begins
                start_time = time.time()

                with open("input_stdin_pyhmmer.fa", "w", encoding="utf-8") as handle:
                    handle.write(target_sequence)
            
                hmm_file = plan7.HMMFile(query_hmmfile)
                seq_file = easel.SequenceFile("input_stdin_pyhmmer.fa", digital=True, alphabet=easel.Alphabet.amino())
                hits_list = list(hmmer.hmmsearch(hmm_file, seq_file, cpus=80))

                pyhmmer_search_fp = open("result_stdin_pyhmmer.domtab", "wb")
                hits_list[0].write(pyhmmer_search_fp, format="domains", header=True)
                for hit in hits_list[1:]:
                    hit.write(pyhmmer_search_fp, format="domains", header=False)
                pyhmmer_search_fp.close()

                search_time = (time.time() - start_time)
                pyhmmsearch_time_2 += search_time
                logging.info(f"pyhmmsearch_time_2: {pyhmmsearch_time_2}")
                # pyhmmer search ends

        except OSError:
            return []
        if not run_result.successful():
            logging.error('hmmsearch returned %d: %s while searching %s',
                          run_result.return_code, run_result.stderr, query_hmmfile)
            raise RuntimeError("Running hmmsearch failed.")
        return list(SearchIO.parse("result.domtab", 'hmmsearch3-domtab'))


def run_hmmsearch_version() -> str:
    """ Get the version of the hmmsearch """

    hmmsearch = get_config().executables.hmmsearch
    command = [
        hmmsearch,
        "-h",
    ]

    help_text = execute(command).stdout
    if not help_text.startswith("# hmmsearch"):
        msg = "unexpected output from hmmsearch: %s, check path"
        raise RuntimeError(msg % hmmsearch)

    version_line = help_text.split('\n')[1]
    return version_line.split()[2]
