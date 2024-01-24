import re

# Function to extract feature information based on conversion type
def extract_feature_info(featureConversionType, line):
    if featureConversionType == 't.to.tName':
        # Define the regex patterns if matching transcript id to transcript name
        feature_id_pattern = r'transcript_id\s+"([^"]+)"'
        feature_name_pattern = r'transcript_name\s+"([^"]+)"'
    elif featureConversionType == 'g.to.gName':
        # Define the regex patterns if matching gene id to gene name
        feature_id_pattern = r'gene_id\s+"([^"]+)"'
        feature_name_pattern = r'gene_name\s+"([^"]+)"'
    elif featureConversionType == 't.to.gName':
        # Define the regex patterns if matching transcript id to gene name
        feature_id_pattern = r'transcript_id\s+"([^"]+)"'
        feature_name_pattern = r'gene_name\s+"([^"]+)"'

    # Search for feature ID and feature name within the line
    feature_id_match = re.search(feature_id_pattern, line)
    feature_name_match = re.search(feature_name_pattern, line)

    if feature_id_match and feature_name_match:
        # Extract feature ID and feature name from the match
        feature_id = feature_id_match.group(1)
        feature_name = feature_name_match.group(1)
        return feature_id, feature_name
    else:
        return None, None

# Function to build a hashmap of feature symbol and name key-value pairs
def build_feature_map(featureConversionType, gtf_file):
    feature_map = {}  # Initialize the hashmap to store feature ID-name key-value pairs

    with open(gtf_file, 'r') as file:
        for line in file:
            # Extract transcript ID and transcript name from each line
            feature_id, feature_name = extract_feature_info(featureConversionType, line)
            if feature_id and feature_name:
                # Add gene ID and gene name as key-value pair to the hashmap
                feature_map[feature_id] = feature_name

    return feature_map

# Function that checks whether the feature_id argument is a novel isoform, then generates appropriate output
def CheckNovelIsoform(feature_id, ConvType, gtfFilePath, novelIsoTitle):
    # This if statement checks whether the conversion type is transcript id to gene symbol, 
    # in which case it returns the feature_name containing the gene symbol (this should happen regardless of novelIsoTitle value)
    if '-' in feature_id and ConvType == 't.to.gName':
        novel_isoform_id = feature_id
        # Change the feature_id to gene id by extracting the substring before the first hyphen (for generating hashmap)
        feature_id = feature_id.split('-')[0]  
        # Build hashmap of gene symbol and name pairs and find the novel isoform's corresponding gene name
        feature_map = build_feature_map('g.to.gName', gtfFilePath)
        feature_name = feature_map.get(feature_id)
        
        if feature_name:
            print(f"Transcript ID: {novel_isoform_id}, Gene Symbol: {feature_name}")
            return feature_name
        else:
            print(f"No matching genes found for novel Isoform ID: {feature_id}")
            return novel_isoform_id
    
    # This elif statement will return the novel isoform id as it is not a gene id, its a transcript id
    # It should return the same thing regardless of novelIsoTitle
    elif '-' in feature_id and ConvType == 'g.to.gName':
        print(f"{feature_id} is not a gene, did you mean ConvType = 't.to.gName?")
        return feature_id
    
    # Then novelIsoformTitle is set to true, check if feature_id is a novel isoform and return appropriate strings
    elif '-' in feature_id and ConvType == 't.to.tName':
        
        # If you don't want to generate the full isoform title, just return the full novel transcript id
        if novelIsoTitle == False:
            print(f"Novel Isoform detected: {feature_id}")
            return feature_id
        
        # The following code in this elif statement assumes novelIsoTitle is set to True
        novel_isoform_id = feature_id
        # Change the feature_id to gene id by extracting the substring before the first hyphen (for generating hashmap)
        feature_id = feature_id.split('-')[0]  
        # Build hashmap of gene symbol and name pairs and find the novel isoform's corresponding gene name
        feature_map = build_feature_map('g.to.gName', gtfFilePath)
        feature_name = feature_map.get(feature_id)
        
        # If novelIsoTitle argument is set to true, return the full Novel Isoform title
        if feature_name:
            fullIsoTitle = f"Novel {feature_name} Isoform \n({novel_isoform_id})"
            print(fullIsoTitle)
            return fullIsoTitle
    
    # else return a boolean false value
    else:
        #print("Feature ID was not detected as a novel isoform")
        return 'No novel isoform detected'
        
# Main function that retrieve the desired corresponding feature name to the provided gene/transcript ID
def FeatureIDtoName(feature_id, ConvType, gtfFilePath, novelIsoTitle = False):
    # Check if conversion type parameter is in correct format
    valid_ConversionTypes = ['t.to.tName', 'g.to.gName', 't.to.gName']
    if ConvType not in valid_ConversionTypes:
        print("Error: Invalid Conversion Type \nRequires: 't.to.tName', 'g.to.gName' or 't.to.gName'")
        exit()
    
    novel_isoform_output = CheckNovelIsoform(feature_id, ConvType, gtfFilePath, novelIsoTitle)
    # Check if feature_id is a novel mrna isoform, and return proper outputs
    if novel_isoform_output != 'No novel isoform detected' :
        return novel_isoform_output

    # This else statement is run if the feature_id wasn't detected as a novel isoform
    else:
        # Build a hashmap of feature symbol and name key-value pairs
        feature_map = build_feature_map(ConvType, gtfFilePath)
        # Retrieve the feature name based on feature ID
        feature_name = feature_map.get(feature_id)
    
        if feature_name:
            if ConvType == 't.to.tName':
                print(f"Transcript ID: {feature_id}, Transcript Symbol: {feature_name}")
            elif ConvType == 'g.to.gName':
                print(f"Gene ID: {feature_id}, Gene Symbol: {feature_name}")
            elif ConvType == 't.to.gName':
                print(f"Transcript ID: {feature_id}, Gene Symbol: {feature_name}")
            return feature_name
        else:
            print(f"No matching names found for the provided feature ID: {feature_id}")
            return feature_id

# Function to get a list of gene names for a given list of gene IDs
def ListofGeneIDstoNames(gene_ids, gtfFilePath):
    gene_names = []  # Initialize the list to store gene names

    # Build a hashmap of gene ID to gene name key-value pairs
    feature_map = build_feature_map('g.to.gName', gtfFilePath)

    for gene_id in gene_ids:
        # Remove the "-PAR-Y" suffix using regular expression
        gene_id_without_suffix = re.sub(r"-PAR-Y$", "", gene_id)
        
        # Retrieve the gene name from the hashmap based on gene ID
        gene_name = feature_map.get(gene_id_without_suffix)
        if gene_name:
            gene_names.append(gene_name)
        else:
            print(f"Couldn't find gene name for {gene_id}, appending None")
            gene_names.append(None)  # Append None if no matching gene name is found

    return gene_names

# Function to get a list of gene IDs for a given list of gene names
def ListofGeneNamestoIDs(gene_names, gtfFilePath):
    gene_ids = []  # Initialize the list to store gene IDs

    # Build a hashmap of gene name to gene ID key-value pairs
    feature_map = build_feature_map('g.to.gName', gtfFilePath)

    for gene_name in gene_names:
        gene_id = None
        # Search the feature_map for the gene ID that matches the provided gene name
        for feature_id, feature_name in feature_map.items():
            if feature_name == gene_name:
                gene_id = feature_id
                break

        if gene_id:
            gene_ids.append(gene_id)
        else:
            print(f"Couldn't find gene ID for gene name: {gene_name}, appending None")
            gene_ids.append(None)  # Append None if no matching gene ID is found

    return gene_ids

def ListofTranscriptIDstoNames(transcript_ids, gtfFilePath):
    transcript_names = []  # Initialize the list to store transcript names

    # Build a hashmap of transcript ID to transcript name key-value pairs
    feature_map = build_feature_map('t.to.tName', gtfFilePath)

    for transcript_id in transcript_ids:
        # Retrieve the transcript name from the hashmap based on transcript ID
        transcript_name = feature_map.get(transcript_id)
        if transcript_name:
            transcript_names.append(transcript_name)
        else:
            print(f"Couldn't find transcript name for {transcript_id}, appending None")
            transcript_names.append(None)  # Append None if no matching transcript name is found

    return transcript_names


def TestScript():
    gtf_file = '/data/gpfs/projects/punim0646/manveer/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
    gene_id = 'ENSG00000104435.14'   # Gene Name: STMN2
    gene_id2 = 'ENSG00000237613.2'   # Gene Name: FAM138A
    gene_list = [gene_id, gene_id2]
    transcript_id = 'ENSG00000104435.14-79611117-79665011-1'
    transcript_id2 = 'ENST00000641515.2' # transcript name : OR4F5-201
    
    geneNameList = ['STMN2', 'FAM138A']
    gtf_file = '/data/gpfs/projects/punim0646/manveer/gencode.v43.chr_patch_hapl_scaff.annotation.gtf'
    print(ListofGeneNamestoIDs(geneNameList, gtf_file))
    
    print("TESTING GENE ID LIST TO NAMES FUNCTION")
    print(ListofGeneIDstoNames(gene_list, gtf_file))
    print('\n\n')
    # Retrieve desired name based on feature symbol and conversion type
    print(f"return value: {FeatureIDtoName(gene_id, 'g.to.gName', gtf_file, novelIsoTitle=True)} \nExpected output: STMN2\n")
    print(f"return value: {FeatureIDtoName(transcript_id, 'g.to.gName', gtf_file, novelIsoTitle=True)} \nExpected output: novel trans ID\n")
    print(f"return value: {FeatureIDtoName(transcript_id, 't.to.tName', gtf_file, novelIsoTitle=True)} \nExpected output: novel trans ID title\n")
    print(f"return value: {FeatureIDtoName(transcript_id, 't.to.gName', gtf_file, novelIsoTitle=True)} \nExpected output: STMN2\n")

    print("Transcript 2 : OR4F5-201--------")
    print(f"return value: {FeatureIDtoName(transcript_id2, 'g.to.gName', gtf_file, novelIsoTitle=True)} \nExpected output: transcript2 ID\n")
    print(f"return value: {FeatureIDtoName(transcript_id2, 't.to.tName', gtf_file, novelIsoTitle=True)} \nExpected output: OR4F5-201\n")
    print(f"return value: {FeatureIDtoName(transcript_id2, 't.to.gName', gtf_file, novelIsoTitle=True)} \nExpected output: OR4F5\n")

    print('\nnovelIsoTitle = False')
    print(f"return value: {FeatureIDtoName(gene_id, 'g.to.gName', gtf_file)} \nExpected output: STMN2\n")
    print(f"return value: {FeatureIDtoName(transcript_id, 'g.to.gName', gtf_file)} \nExpected output: Novel trans ID\n")
    print(f"return value: {FeatureIDtoName(transcript_id, 't.to.tName', gtf_file)} \nExpected output: Novel trans ID\n")
    print(f"return value: {FeatureIDtoName(transcript_id, 't.to.gName', gtf_file)} \nExpected output: STMN2\n")

    print("Transcript 2 : OR4F5-201--------")
    print(f"return value: {FeatureIDtoName(transcript_id2, 'g.to.gName', gtf_file)} \nExpected output: transcript2 ID\n")
    print(f"return value: {FeatureIDtoName(transcript_id2, 't.to.tName', gtf_file)} \nExpected output: OR4F5-201\n")
    print(f"return value: {FeatureIDtoName(transcript_id2, 't.to.gName', gtf_file)} \nExpected output: OR4F5\n")