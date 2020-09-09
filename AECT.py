import random
import numpy as np
import pandas as pd
import tensorflow as tf
import argparse

from ae.preprocess import load_file, preprocess
from ae.autoencoder import Autoencoder, train

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        usage="%prog -i <input_data> -f [file_type] -o [output_data]",
        description='Autoencoder on Cell-free DNA TSS coverage profile',
        epilog = "Written by Han Bowei (hanbw0120@foxmail.com), 2020\n"
    )

    parser.add_argument('--input_data', '-i', type=str, help='Input dataset path', required=True)
    parser.add_argument('--output_data', '-o', type=str, default=None, help='Output dataset path (default:input_data.imputed.csv)')
    parser.add_argument('--file_type', '-f', type=str, choices=['csv', 'tab'], default='csv',
                        help='Input dataset type (csv or tab), (default: csv)')
    parser.add_argument('--seed', type=int, default=0, help='Random seed for reproducing')
    parser.add_argument('--threads', '-t', type=int, default=None, help='Number of threads for training (default: all cores)')

    # data preprocess options
    parser.add_argument('--log1p', '-l', type=str, choices=['True', 'False'], default='True',
                        help='Log transformed data (default: True)')
    parser.add_argument('--scale', '-s', type=str, choices=['True', 'False'], default='True',
                        help='Scale data (default: True)')

    # network and training options
    parser.add_argument('--hidden_size', type=str,  default='128,64,32,64,128',
                        help='Hidden layers, separated with commas (default: 128,64,32,64,128)')
    parser.add_argument('--batch_size', '-b', type=int, default=32, help='Batch size (default: 32)')
    parser.add_argument('--epochs', '-e', type=int, default=500, help='Epochs (default: 500)')
    parser.add_argument('--reduce_lr', type=int, default=10, help='reduce_lr (default: 10)')
    parser.add_argument('--early_stop', type=int, default=15, help='early_stop (default: 15)')

    args = parser.parse_args()

    # Set random seed
    seed = args.seed
    random.seed(seed)
    np.random.seed(seed)
    tf.set_random_seed(seed)
    tf.compat.v1.set_random_seed(seed)
    
    # Load data
    dataset = args.input_data
    file_type = args.file_type
    data_in = load_file(dataset, file_type = file_type)

    # Pre-process data
    log1p_data = True if args.log1p == 'True' else False
    scale_data = True if args.scale == 'True' else False
    data_process = preprocess(data_in, log1p=log1p_data, scale_data= scale_data)

    # Build networks
    hidden_size = args.hidden_size.strip().split(',')
    hidden_size = list(map(int, hidden_size))
    network_kwds = {'input_size': data_in.shape[1], 'hidden_size': hidden_size}
    net = Autoencoder(**network_kwds)
    net.build()

    # Training
    train(data_process, data_in, net, epochs=args.epochs, reduce_lr=args.reduce_lr,
        early_stop=args.early_stop, batch_size=args.batch_size, threads=args.threads)

    # Prediction
    data_predict = net.predict(data_process)

    # Generating Output
    data_out = pd.DataFrame(data_predict)
    data_out.columns = data_in.columns
    data_out.index = data_in.index

    if args.output_data:
        output_file_name = args.output_data
    elif file_type == 'csv':
        output_file_name = dataset.replace('.csv', '') + '.imputed.csv'
    else:
        output_file_name = dataset + '.imputed.csv'
    pd.DataFrame(data_out).T.to_csv(output_file_name, float_format='%.3f')
