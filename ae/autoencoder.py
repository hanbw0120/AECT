import numpy as np
import tensorflow as tf
import keras
from keras.layers import Input, Dense, Activation, BatchNormalization
from keras import backend as K
from keras.models import Model


class Autoencoder():
    def __init__(self, input_size, output_size=None, hidden_size=[128, 64, 32, 64, 128]):
        self.input_size = input_size
        if output_size:
            self.output_size = output_size
        else:
            self.output_size = input_size
        self.hidden_size = hidden_size
        self.input_layer = None
        self.model = None

    def build(self):
        self.input_layer = Input(shape=(self.input_size,))
        last_hidden = self.input_layer

        for neuron_count in self.hidden_size:
            last_hidden = Dense(neuron_count, activation=None, kernel_initializer='glorot_uniform')(last_hidden)
            last_hidden = BatchNormalization(center=True, scale=False)(last_hidden)
            last_hidden = Activation('relu')(last_hidden)

        output_layer = Dense(self.output_size, kernel_initializer='glorot_uniform')(last_hidden)

        self.model = Model(inputs=self.input_layer, outputs=output_layer)

    def predict(self, data_process):
        print('Generating imputed data...')
        data_predict = self.model.predict(data_process)
        return data_predict


def train(inputs, outputs, network, epochs=500, reduce_lr=10, early_stop=15,
          batch_size=32, clip_grad=5., validation_split=0.2, threads=None):
    K.set_session(tf.Session(config=tf.ConfigProto(intra_op_parallelism_threads=threads, inter_op_parallelism_threads=threads)))
    model = network.model
    loss = keras.losses.mean_squared_error
    optimizer = keras.optimizers.RMSprop(clipvalue=clip_grad)
    model.compile(loss=loss, optimizer=optimizer)

    callbacks = []
    if reduce_lr:
        reduce_lr_callback = keras.callbacks.ReduceLROnPlateau(monitor='val_loss', patience=reduce_lr, verbose = False)
        callbacks.append(reduce_lr_callback)
    if early_stop:
        early_stop_callback = keras.callbacks.EarlyStopping(monitor='val_loss', patience=early_stop, verbose = False)
        callbacks.append(early_stop_callback)

    model.fit(inputs, outputs, epochs=epochs, batch_size=batch_size,
              shuffle=True, callbacks=callbacks, validation_split=validation_split)
