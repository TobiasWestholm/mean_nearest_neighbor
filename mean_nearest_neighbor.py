# coding: utf-8
# author: Tobias Westholm
# email: tobias_westholm@hotmail.com
import pandas as pd

def mean_nearest_distance(df, image, ROI):
    """
    simple function that calculates and returns the mean shortest distance 
    for cell type 2 to cell type 1 for the specified ROI ID.
    """
    if df.empty:
        # cannot calculate distance if there are no cells of type 2
        mean_distance_to_type1 = "NA"
    else:
        # calculate mean nearest distance from cell type 2 to cell type 1
        mean_distance_to_type1 = df["distance_to_type1"].mean()

    # create results list
    results = [image, ROI, mean_distance_to_type1]

    return results
  
if __name__ == "__main__":
    #VARIABLES
    in_path = "C:\\Users\\tobia\\Desktop\\Exjobb\\Python_stuff\\nearest_neighbor\\ROI_centroid_distances.tsv"       #path to input csv file
    out_path = "C:\\Users\\tobia\\Desktop\\Exjobb\\Python_stuff\\nearest_neighbor"      #path to the results directory, default is current directory'
    pair = "CD45:PANCK" #cell type pair
    separator = '\t'    #csv separator, default tab
    decimal = '.'       #float decimal sign, default: .
    
    # create the output dictionary
    out = {'image':[],
    'ROI':[],
    'type2_nearest_celltype1':[]}

    # read the csv file
    df = pd.read_csv(in_path,
                     sep=separator,
                     decimal=decimal,
                     low_memory=False,
                     skiprows=1,
                     usecols=[0, 1, 2, 3],
                     names=['Image',
                            'Class',
                            'ROI',
                            'distance_to_type1']).dropna(axis=1)
    #define cells
    cell_type2 = pair.split(':')[1]
    cell_type1 = pair.split(':')[0]
    # list all images
    images = df['Image'].unique()
    # check that all images have both classes
    for image in images:
        if cell_type1 not in df.groupby('Image').get_group(
            image)['Class'].unique():
            print(f"{cell_type1} not in {image}")
        if cell_type2 not in df.groupby('Image').get_group(
            image)['Class'].unique():
            print(f"{cell_type2} not in {image}")
    # iterate over images
    for image in images:
        # filter on cell type included in analysis AND clean out any annotations & detections (outside of ROIs) with the parent named "Image"
        filtered = df[(df['Class'] == cell_type1) | (
                  df['Class'] == cell_type2) & (df['ROI'] != 'Image')].reset_index(drop=True)
        # extract the image
        one_pic = filtered.groupby('Image').get_group(image).reset_index(drop=True)
        # iterate through all ROIs in the picture
        for roi_value in sorted(one_pic['ROI'].unique()):
            # pick out all rows with the right ROI number. roi_value MIGHT NEED TO BE CONVERTED TO STRING OR INTEGER IN THE COMPARISON BELOW
            one_roi = one_pic.loc[one_pic['ROI'] == roi_value].reset_index(drop=True)
            one_roi = one_roi.drop(one_roi[one_roi.Class == cell_type1].index)
            try:
            # calculate mean nearest distances for each cell of type 2 to type 1.
                results = mean_nearest_distance(one_roi, image, roi_value)
            except Exception as e:
                print(f'Error in {image}: {e}')
                continue
            # append to output dictionary
            for key, value in zip(out.keys(), results):
                out[key].append(value)
    # create dataframe from out dictionary
    out_df = pd.DataFrame(out)
    print(out_df.head())

    # save the output
    out_df.to_csv(
        f'{out_path}/{cell_type2}_mean_nearest_distance_to_{cell_type1}.csv',
        index=False,
        sep='\t')
