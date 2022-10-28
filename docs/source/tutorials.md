Transcriptome reference generation
----------------------------------
.. _ref_gen:

The first steps of cerberus are to generate a reference set of transcript start
sites (TSSs), intron chains (ICs), and transcript end sites (TESs). TSSs and TESs
can either be added from GTFs or BED files. ICs can be added from GTFs. To extract
triplet features from a GTF, see :ref:`this <gtf_feats>` section.


Extracting triplet features from GTFs
-------------------------------------
.. _gtf_feats:

To extract TSSs or TESs from a GTF and save them in BED file format, use ``cerberus gtf_to_bed``.
