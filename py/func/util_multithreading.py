from joblib import Parallel, delayed
import collections
import pandas as pd
import numpy as np
import segmentImporter



# def helper_map(df,conv_dict,seg):
#     time.sleep(.5)
#     out = df[seg].map(lambda x: test_conversion(x,conv_dict))
#     return out

# def applyParallel(df, func, ncore, conv_dict, seg):
#     retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(func)(df_chunk,conv_dict,seg) for df_chunk in df)
#     return pd.concat(retLst)



####################################
# Read segmentation
# Multithreading implementation for read segmentation
# def segmentation_parallel_wrapper(fastq_chunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set,ncore):
#     parsedSeqDict=collections.defaultdict(list)
#     parsedQualDict=collections.defaultdict(list)
#     Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(
#         delayed(segmentImporter.fastq_segmentation)(parsedSeqDict,parsedQualDict,subchunk,settings,headerSplitRegex,readKey,segment_parsed,segment_parsed_set) 
#             for subchunk in segmentImporter.split_fastq_per_lines(fastq_chunk,n=40000)) # Each CPU processes 40000 FASTQ records
#     return parsedSeqDict,parsedQualDict



# ####################################
# # Quality filtering
# # Multithreading implementation for quality filtering
# def qfilter_parallel_wrapper(seq_chunk, qual_chunk, qc_targets, qscore_dict, ncore):
#     # Further split the data chunks into subchunks by the number of CPUs
#     seq_subchunks = np.array_split(seq_chunk,ncore)
#     qual_subchunks = np.array_split(qual_chunk,ncore)

#     retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(segmentImporter.qualityFilteringForDataFrame)(df_set, qc_targets, qscore_dict) for df_set in zip(seq_subchunks,qual_subchunks))
#     return pd.concat(retLst)



####################################
# Error correction
# # Helper function for error correction
# def sequenceCorrection(seq_chunk,suggestion_verbosity,correctOpt,correction_method,symspelldb,wlset):
#     seq_chunk = seq_chunk.map(lambda x:findClosestCandidate(x,suggestion_verbosity,correctOpt[correction_method]["levenshtein_distance"],symspelldb,wlset))
#     return seq_chunk

# # Return closest hit
# def findClosestCandidate(seq,suggestion_verbosity,dist_thresh,symspelldb,wlset):
#     if wlset:
#         if seq in wlset:
#             return seq
#     suggestions = symspelldb.lookup(seq,suggestion_verbosity,dist_thresh,transfer_casing=True)
#     if len(suggestions)>=1:
#         top_candidate=suggestions[0].term.upper()
#         return top_candidate
#     else:
#         return "-"

# # Multithreading implementation for error correction
# def correct_parallel_wrapper(seq_series,suggestion_verbosity,correctOpt,correction_method,symspelldb,ncore,wlset=None):
#     # Split the sequence series into chunks by the number of CPUs
#     seq_chunks = np.array_split(seq_series,ncore)
    
#     retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(sequenceCorrection)(seq_chunk,suggestion_verbosity,correctOpt,correction_method,symspelldb,wlset) for seq_chunk in seq_chunks)
#     return pd.concat(retLst)



####################################
# Create a clean segment sequence table / Value table generation
# # Create a clean segment sequence table
# def gen_clean_segment_seq(df_subchunk,segments_raw,correctionDictionaries,settings):
#     cols_ordered=["Header"]+[i for i in sorted(df_subchunk.columns) if not i=="Header"]
#     cols_ordered_final=[]

#     for ncol,segment_raw_now in enumerate(cols_ordered):
#         if ncol>0:
#             if not segment_raw_now in segments_raw:
#                 df_subchunk=df_subchunk.drop(segment_raw_now,axis=1)
#                 continue
#             col_index=segments_raw.index(segment_raw_now)
#             component_corrected_now=settings.corrected_components[col_index]
#             opt_now=settings.correctOptDict[component_corrected_now]
            
#             if (component_corrected_now in correctionDictionaries) and ("I2M_CORRECTION" in opt_now["func_ordered"] or "M2A_CORRECTION" in opt_now["func_ordered"] or "BARTENDER_CORRECTION" in opt_now["func_ordered"] or "CUSTOM_CORRECTION" in opt_now["func_ordered"]):
#                 df_subchunk[segment_raw_now]=df_subchunk[segment_raw_now].map(lambda x: seq_correct_and_write(x,reference=correctionDictionaries[component_corrected_now]["correctionDict"]))
#             df_subchunk.rename(columns={segment_raw_now:segment_raw_now+":"+component_corrected_now},inplace=True)
#             cols_ordered_final.append(segment_raw_now+":"+component_corrected_now)
#         else:
#             cols_ordered_final.append("Header")

#     df_subchunk = df_subchunk[cols_ordered_final]
#     return df_subchunk

# # Multithreading implementation for cleas segment sequence table generation
# def qfilter_parallel_wrapper(df_chunk, qual_chunk, qc_targets, qscore_dict, ncore):
#     # Further split the data chunks into subchunks by the number of CPUs
#     df_subchunks = np.array_split(df_chunk,ncore)

#     retLst = Parallel(n_jobs=ncore,require='sharedmem',verbose=10)(delayed(segmentImporter.qualityFilteringForDataFrame)(df_set, qc_targets, qscore_dict) for df_set in zip(seq_subchunks,qual_subchunks))
#     return pd.concat(retLst)
