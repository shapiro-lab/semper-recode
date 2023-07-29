=====
Usage
=====

.. .. include:: ../semper_recode/semper_recode.py

.. .. automodule:: semper_recode
    :members:

Overview:
----------------

- :func:`process_sequence`
    - :func:`modify_TIS_in_frame`
        - :func:`find_in_frame`
        - :func:`efficiency_level`
    - :func:`modify_TIS_out_of_frame`
        - :func:`find_out_of_frame_list`
        - :func:`find_out_of_frame_index`
        - :func:`return_key`
        - :func:`get_alternative_codon`

Archived
    - :func:`to_fasta`
    - :func:`filtered_sequence_eff`
    - :func:`find_lower_eff_sequence`

------------

Global Variables
----------------

    .. py:data:: START_CODON

        A list of canonical start codons. Future versions may include non-canonical start codons and will be appended to this list.

    .. py:data:: MASTER_DF

        The master dataframe containing data required for the analysis such as TIS sequence, efficiency level, and amino acids associated with the sequence.

    .. py:data:: CODON_LIST

        The codon list loaded from a .pkl file.

    .. py:data:: CODON_DICT

        The codon dictionary loaded from a .pkl file.

    .. py:data:: FOUR_LETTERS_CODE

        List of unique four-letter codes from the master dataframe.

------------

Functions
---------

.. py:method:: process_sequence(self)
    
        Return the modified sequence with lower efficiency (if possible) of the nucleotide sequence.
        Loop through modify_TIS_in_frame() and modify_TIS_out_of_frame() to ensure there's no 
        out-of-frame internal AUG and that the expressional level is the lowest possible

    :returns: The modified sequence with lower TIS efficiency (if any).
    :rtype: str

------------

.. py:method:: modify_TIS_in_frame(self, sequence)
    
        Takes in sequence and return modified sequence with lower TIS efficiency by getting the
        location of internal AUG (index) from find_in_frame() and modify the TIS -6 to + 3 of each index
        by replacing the TIS sequence with the sequence with lowest possible efficiency level according 
        to master dataframe
        
   :param sequence: The input nucleotide sequence.
   :type sequence: str

   :returns: The modified sequence with lower TIS efficiency.
   :rtype: str

------------

.. py:method:: find_in_frame(self, sequence)
    
        Takes in a sequence (str) and returns a list of indices where the AUG codon is found.

   :param sequence: The input nucleotide sequence.
   :type sequence: str

   :returns: A list of integers representing the indices of in-frame AUG codons in the sequence.
   :rtype: list

------------

.. py:method:: efficiency_level(self, sequence)
    
        Takes in sequence and return efficiency level

   :param sequence: The input nucleotide sequence.
   :type sequence: str

   :returns: The efficiency level of the sequence.
   :rtype: int

------------

.. py:method:: modify_TIS_out_of_frame(self, sequence)
    
        Takes in sequence (str) and return modified sequence
        with lower TIS efficiency (str)

   :param sequence: The input nucleotide sequence.
   :type sequence: str

   :returns: The modified sequence with lower TIS efficiency.
   :rtype: str

------------

.. py:method:: return_key(self, codon)
    
        Takes in a codon sequence and returns the abbreviation of the codon.

   :param codon: The input codon sequence.
   :type codon: str

   :returns: The abbreviation of the codon.
   :rtype: str

------------

.. py:method:: find_out_of_frame(self, sequence)
    
        Takes in a sequence and returns a list of indices where the out-of-frame AUG codon is found.

   :param sequence: The input nucleotide sequence.
   :type sequence: str

   :returns: A list of integers representing the indices of out-of-frame AUG codons in the sequence.
   :rtype: list

------------

.. py:method:: get_aa_alternative(self, original_codon)
    
        Takes in an original codon and returns the codon producing the key amino acid with the highest fraction.

   :param original_codon: The original codon sequence.
   :type original_codon: str

   :returns: The codon and the fraction of that codon.
   :rtype: tuple

------------
**Archived**

.. py:method:: to_fasta(self, sequence, output_file_name)

        Takes in a list of modified sequences and converts them back to FASTA format for exporting to users.

   :param sequence: The modified list of sequences to be converted to a FASTA file.
   :type sequence: list

   :param output_file_name: The desired name of the output FASTA file.
   :type output_file_name: str

------------

.. py:method:: filtered_sequence_eff(self, efficiency)

        Takes in the efficiency level input and finds the sequences in the dataframe by filtering,
        then displays all the possible sequences that have an efficiency lower than the given value
        and the percentage of found and missing amino acid sequences in the form of a dataframe.

   :param efficiency: The efficiency level to filter the dataframe.
   :type efficiency: float

   :returns: A dataframe containing the filtered sequences.
   :rtype: pandas.DataFrame

------------

.. py:method:: find_lower_eff_sequence(self, efficiency, seq)

        Takes in a nucleotide sequence and checks if there is another sequence that produces
        the same protein but with a lower efficiency level than the input level.

   :param efficiency: The efficiency level to filter the dataframe.
   :type efficiency: float

   :param seq: The nucleotide sequence.
   :type seq: Bio.Seq

   :raises ValueError: If the given sequence doesn't have replacing sequence with efficiency lower
                       than the input efficiency level.

   :returns: Nucleotide sequence with lower efficiency level.
   :rtype: Bio.Seq

------------

Sample Usage
------------

Here's an example of how to use the SemperRecode class:

.. code-block:: python

   # Import the module
   from semper import SemperRecode

   # Create a SemperRecode object with the user input sequence
   user_sequence = "ATGCATCGATCGATCG"
   semper_obj = SemperRecode(user_sequence)

   # Process the sequence to obtain a modified sequence with lower efficiency
   modified_sequence = semper_obj.process_sequence()

   # Convert the modified sequence to FASTA format and export it to a file
   semper_obj.to_fasta(modified_sequence, "modified_sequence")

   # Find sequences with lower efficiency in the dataframe
   filtered_df = semper_obj.filtered_sequence_eff(0.5)

   # Find a nucleotide sequence with lower efficiency
   lower_eff_sequence = semper_obj.find_lower_eff_sequence(0.4, "ATGCTAGCTAGCTAG")