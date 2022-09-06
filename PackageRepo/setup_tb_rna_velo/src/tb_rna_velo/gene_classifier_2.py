from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from scipy.optimize import fsolve
from math import exp
from pickle import dump, load
import tensorflow as tf
import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import tempfile
import os
import random

# Global params here
mpl.rcParams['figure.figsize'] = (12, 10)
colors = plt.rcParams['axes.prop_cycle'].by_key()['color']


class GeneClassifier_2:
    """Class to classify a gene(s)
    """

    def __init__(self, adata=None, adata2=None):
        self.history = None
        self.test_labels = None
        self.test_features = None
        self.val_features = None
        self.train_features = None
        self.val_labels = None
        self.train_labels = None
        self.adata = adata
        self.class_list = ['bi', 'poo', 'bias', 'tarme', 'z']

        self.preprocess()

        self.METRICS = [
            # tf.keras.metrics.TruePositives(name='tp'),
            # tf.keras.metrics.FalsePositives(name='fp'),
            # tf.keras.metrics.TrueNegatives(name='tn'),
            # tf.keras.metrics.FalseNegatives(name='fn'),
            # tf.keras.metrics.CategoricalAccuracy(),
            tf.keras.metrics.SparseCategoricalAccuracy(),
            # tf.keras.metrics.Precision(name='precision'),
            # tf.keras.metrics.Recall(name='recall'),
            # tf.keras.metrics.AUC(name='auc'),
            # tf.keras.metrics.AUC(name='prc', curve='PR'), # precision-recall curve
        ]

        self.EPOCHS = 1000
        self.BATCH_SIZE = 500

        self.early_stopping = tf.keras.callbacks.EarlyStopping(
            monitor='val_loss',
            verbose=1,
            patience=5,
            restore_best_weights=True)

        self.train_class_counts = {
            'bi': 0,
            'poo': 0,
            'biased': 0,
            'tb_arme': 0,
            'zero_inflated': 0,
            'total': 0}

        self.train_dataset = None
        self.test_dataset = None
        self.bias_init = None
        self.scaler_path = "data/scaler_2.pkl"
        self.scaler = None
        self.model = None
        self.model_path = "data/model/model_2.h5"
        self.checkpoint_path = "data/checkpoints/cp_2.ckpt"
        # self.label_data()
        # self.get_aggregate_value()

    def preprocess(self):
        pass

    def label_data(self):
        """Method to label data into required classes based on p_values and standard deviation
        """

        # zero_inflated = self.adata.var[(self.adata.var.sum_allele_1 <= 50) & (self.adata.var.sum_allele_2 <= 50)]
        # zero_inflated_index = zero_inflated.index
        # print(len(zero_inflated_index))

        # non_zero_inflated = self.adata.var[(self.adata.var.sum_allele_1 > 50) & (self.adata.var.sum_allele_2 > 50)]


        # use adata p_value, ratio and std values to give a particular class label to each gene
        group_1 = self.adata.var[self.adata.var.p_value >= 0.00000005]
        group_2 = self.adata.var[self.adata.var.p_value < 0.00000005]

        # group_1 = non_zero_inflated[non_zero_inflated.p_value >= 0.00000005]
        # group_2 = non_zero_inflated[non_zero_inflated.p_value < 0.00000005]

        bi_allelic = group_1[((group_1.ratio_std_allele_1 >= 0.1) & (group_1.ratio_std_allele_1 <= 0.4)) &
                             ((group_1.sum_allele_1 + group_1.sum_allele_2) > 80)
                             ]
        bi_allelic_index = bi_allelic.index
        print(len(bi_allelic_index))

        parent_of_origin = group_2[((group_2.ratio_allele_1 <= 0.1) |
                                   (group_2.ratio_allele_1 >= 0.9) |
                                   (group_2.ratio_allele_2 <= 0.1) |
                                   (group_2.ratio_allele_2 >= 0.9)) &
                                   ((group_2.sum_allele_1 + group_2.sum_allele_2) > 80)
                                   ]
        parent_of_origin_index = parent_of_origin.index
        print(len(parent_of_origin_index))

        other_1_index = set([x for x in list(group_1.index) if x not in list(bi_allelic_index)])
        other_2_index = set([x for x in list(group_2.index) if x not in list(parent_of_origin)])
        other_index = other_1_index.union(other_2_index)
        other = self.adata.var[self.adata.var.index.isin(other_index)]

        biased = other[
            ((other.p_value < 0.00000005) & (other.ratio_std_allele_1 <= 0.4)) &
            ((other.sum_allele_1 + other.sum_allele_2) > 80)
            ]
        biased_index = biased.index
        print(len(biased_index))

        other_index = set([x for x in list(other_index) if x not in list(biased_index)])
        other = self.adata.var[self.adata.var.index.isin(other_index)]

        tb_arme_1 = other[(other.p_value < 0.00000005) &
                          (
                                  ((other.ratio_allele_1 > 0.1) & (other.ratio_allele_1 < 0.9)) |
                                  ((other.ratio_allele_2 > 0.1) & (other.ratio_allele_2 < 0.9))
                          ) &
                          ((other.sum_allele_1 + other.sum_allele_2) > 80)
                          ]
        tb_arme_2 = other[(other.p_value >= 0.00000005) &
                          ((other.ratio_std_allele_1 > 0.4) & (other.ratio_std_allele_1 < 0.6)) &
                          ((other.sum_allele_1 + other.sum_allele_2) > 80)
                          ]
        tb_arme = pd.concat([tb_arme_1, tb_arme_2])
        tb_arme_index = tb_arme.index
        print(len(tb_arme_index))

        zero_inflated_index = [x for x in list(other_index) if x not in list(tb_arme_index)]
        zero_inflated = self.adata.var[self.adata.var.index.isin(zero_inflated_index)]
        print(len(zero_inflated_index))


        label = []
        class_name = []
        for x in self.adata.var.index:
            if x in bi_allelic_index:
                class_name.append('bi')
                label.append(0)
            elif x in parent_of_origin_index:
                class_name.append('poo')
                label.append(1)
                
            elif x in biased_index:
                class_name.append('biased')
                label.append(2)
                
            elif x in tb_arme_index:
                class_name.append('tarme')
                label.append(3)
                
            elif x in zero_inflated_index:
                class_name.append('z')
                label.append(4)
        #     else:
        #         label.append(2)
        self.adata.var['class_names'] = class_name
        self.adata.var['class_label'] = label
        
        return

    def get_aggregate_value(self):
        """Method to calculate the aggregate values for the count matrices
        """
        allele_1 = pd.DataFrame(self.adata.layers['spliced'], columns=self.adata.var.index)
        allele_2 = pd.DataFrame(self.adata.layers['unspliced'], columns=self.adata.var.index)
        allele_1_T = allele_1.transpose(copy=True)
        allele_2_T = allele_2.transpose(copy=True)
        allele_T = (allele_1_T - allele_2_T) / (allele_1_T + allele_2_T)
        self.adata.layers['aggregate'] = allele_T.transpose(copy=True)
        return

    def generate_train_test_set(self):
        """Method to generate training and test dataset using the aggregated value
        """
        aggregate_layer = self.adata.layers['aggregate']
        aggregate_layer_T = aggregate_layer.transpose()

        aggregate_df = pd.DataFrame(aggregate_layer_T, index=self.adata.var.index)

        # need to replace NAs with 0
        aggregate_df = aggregate_df.fillna(0)

        aggregate_df['class_names'] = self.adata.var['class_names']
        aggregate_df['class_label'] = self.adata.var['class_label']
        

        bi, poo, biased, tb_arme, zero_inflated = np.bincount(self.adata.var['class_label'])
        print(bi, poo, biased, tb_arme, zero_inflated)

        # TODO- get the test set numbers using the total counts
        test_set_1 = aggregate_df[aggregate_df['class_label'] == 0].head(int(0.2 * bi)).index
        test_set_2 = aggregate_df[aggregate_df['class_label'] == 1].head(int(0.2 * poo)).index
        test_set_3 = aggregate_df[aggregate_df['class_label'] == 2].head(int(0.2 * biased)).index
        test_set_4 = aggregate_df[aggregate_df['class_label'] == 3].head(int(0.2 * tb_arme)).index
        test_set_5 = aggregate_df[aggregate_df['class_label'] == 4].head(int(0.2 * zero_inflated)).index

        # print(test_set_1, test_set_2, test_set_3, test_set_4, test_set_5)
        
        print(test_set_4)

        test_set = test_set_1.union(test_set_2)
        test_set = test_set.union(test_set_3)
        test_set = test_set.union(test_set_4)
        test_set = test_set.union(test_set_5)
        test_set = list(set(test_set))

        test_set_shuffled = random.sample(test_set, len(test_set))

        train_set = [x for x in list(self.adata.var.index) if x not in test_set]
        train_set_shuffled = random.sample(train_set, len(train_set))

        # TODO: option to save test sets as .csv file if the path is given

        self.train_dataset = aggregate_df.loc[train_set_shuffled]

        self.test_dataset = aggregate_df.loc[test_set_shuffled]
        # shuffle columns. TODO: find better way to shuffle
#         label_column = self.test_dataset['class_label']
#         class_name_column = self.test_dataset['class_name_label']
        
#         self.test_dataset = self.test_dataset.drop('class_label', axis=1)
#         self.test_dataset = (self.test_dataset.transpose().sample(frac=1)).transpose()
#         self.test_dataset['class_label'] = label_column

        return

    def get_test_dataset(self):
        return self.test_dataset

    # TODO: Combine generate_train_features and generate_test_features
    def generate_train_features(self, load_scaler: bool = False):
        """generate training features for the model to process

        Args:
            load_scaler (bool, optional): parameter to load already saved Scaler model. Defaults to False.
        """
        # separate into train and validation sets
        # self.train_dataset = self.train_dataset.fillna(0)

        y = self.train_dataset.iloc[:, -1]
        X = self.train_dataset.iloc[:, :-2]
        
        print(self.train_dataset.head(1))

        X_train, X_val, y_train, y_val = train_test_split(X, y, test_size=0.3, random_state=0)

        self.train_features = np.array(X_train)
        self.val_features = np.array(X_val)
        self.train_labels = np.array(y_train)
        self.val_labels = np.array(y_val)
        
        print(set(self.train_labels))

        self.train_class_counts['bi'], self.train_class_counts['poo'], self.train_class_counts['biased'],\
        self.train_class_counts['tb_arme'], self.train_class_counts['zero_inflated'] = np.bincount(self.train_labels)

        self.train_class_counts['total'] = self.train_class_counts['bi'] + self.train_class_counts['poo'] + \
                                           self.train_class_counts['biased'] + self.train_class_counts['tb_arme'] + \
                                           self.train_class_counts['zero_inflated']

        self.get_initial_bias()

        if load_scaler is False:
            self.scaler = StandardScaler()
            self.scaler.fit_transform(self.train_features)
            self.save_scaler(self.scaler_path)
        else:
            self.load_scaler(self.scaler_path)

        self.train_features = self.scaler.transform(self.train_features)
        self.val_features = self.scaler.transform(self.val_features)

        # clipping not necessary as the data would be in the range of -1,1
        # train_features = np.clip(train_features, -5, 5)
        # val_features = np.clip(val_features, -5, 5)

        return

    def generate_test_features(self):
        """generate test features for the model to process
        """
        if self.scaler is None:
            self.load_scaler(self.scaler_path)

        y_test = self.test_dataset.iloc[:, -1]
        X_test = self.test_dataset.iloc[:, :-2]
        self.test_features = np.array(X_test)
        self.test_features = self.scaler.transform(self.test_features)
        # test_features = np.clip(test_features, -5, 5)
        self.test_labels = np.array(y_test)

    def save_scaler(self, scaler_path: str):
        """save scaler object to given path

        Args:
            scaler_path (str): scaler object file path
        """
        dump(self.scaler, open(scaler_path, 'wb'))
        return

    def load_scaler(self, scaler_path: str):
        """load scaler object from given path

        Args:
            scaler_path (str): scaler object file path
        """
        self.scaler = load(open(scaler_path, 'rb'))
        pass

    def load_train_set(self, train_data_path: str):
        """loads training dataset into class from a file (.csv format)

        Args:
            train_data_path (str): training dataset file path
        """
        self.train_dataset = pd.read_csv(train_data_path, sep=";", index_col=[0])
        self.generate_train_features()
        return

    def load_test_set(self, test_data_path: str):
        """loads test dataset into class from a file (.csv format)

        Args:
            test_data_path (str): test dataset file path
        """
        self.test_dataset = pd.read_csv(test_data_path, sep=";", index_col=[0])
        self.generate_test_features()
        return

    def get_initial_bias(self):
        """Generate initial bias to train moder quickly (start with lower loss)
        """
        # define the frequency of different classes
        f = (self.train_class_counts['bi'] / self.train_class_counts['total'],
             self.train_class_counts['poo'] / self.train_class_counts['total'],
             self.train_class_counts['biased'] / self.train_class_counts['total'],
             self.train_class_counts['tb_arme'] / self.train_class_counts['total'],
             self.train_class_counts['zero_inflated'] / self.train_class_counts['total']
             )

        # define the equation
        def eqn(x, frequency):
            sum_exp = sum([exp(x_i) for x_i in x])
            return [exp(x[i]) / sum_exp - frequency[i] for i in range(len(frequency))]

        # calculate bias init
        self.bias_init = list(fsolve(func=eqn, x0=np.array([0] * len(f)), args=(f,)))
        
        print(self.bias_init)


        self.bias_init = tf.keras.initializers.Constant(self.bias_init)
        

        return

    def save_model(self, model_path: str):
        """save the trained model into filesystem

        Args:
            model_path (str): file path for the model (currently only .h5 format)
        """
        tf.keras.models.save_model(self.model, model_path)
        return

    def load_model(self, model_path: str):
        """load the saved model from filesystem

        Args:
            model_path (str): file path for the model (currently only .h5 format)
        """
        self.model = tf.keras.models.load_model(model_path)
        return

#     def make_model(self):
#         """create classification network architecture using fully connected layers
#         TODO: to be combined with or replaced by make_model_conv
#         Returns:
#             tf.keras.Model: classification network model
#         """
#         # design the neural network with two dense layers
#         inputs = tf.keras.layers.Input(shape=(self.train_features.shape[-1],))
#         x = tf.keras.layers.Dense(64, activation="relu", name="dense_1")(inputs)
#         x = tf.keras.layers.Dropout(0.6)(x)
#         x = tf.keras.layers.Dense(32, activation="relu", name="dense_2")(x)
#         x = tf.keras.layers.Dropout(0.6)(x)
#         x = tf.keras.layers.Dense(16, activation="relu", name="dense_3")(x)
#         outputs = tf.keras.layers.Dense(5, name="predictions", bias_initializer=self.bias_init)(x)

#         model = tf.keras.Model(inputs=inputs, outputs=outputs)

#         model.compile(optimizer=tf.keras.optimizers.RMSprop(learning_rate=0.0001),
#                       loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
#                       metrics=self.METRICS)
#         return model

    def make_model_conv(self):
        """create classification network architecture using convolution layers and fully connected layers

        Returns:
            _type_: tf.keras.Model: classification network model
        """
        # design the neural network with two dense layers
        inputs = tf.keras.layers.Input(shape=(self.train_features.shape[-1],))
        
#         x = tf.keras.layers.Dense(128, activation="relu", name="dense_0")(inputs)
#         x = tf.keras.layers.Dropout(0.2)(x)
        

        x = tf.keras.layers.Lambda(lambda y: tf.keras.backend.expand_dims(y, axis=-1))(inputs)
        x = tf.keras.layers.Conv1D(kernel_size=16, filters=64)(x)
        x = tf.keras.layers.MaxPooling1D(pool_size=64)(x)
        x = tf.keras.layers.Dropout(0.25)(x)
        # x = tf.keras.layers.Dropout(0.3)(x)

        x = tf.keras.layers.Flatten()(x)

#         # non-linear layer for attention
#         # TODO: calculate the attention from paper
        x = tf.keras.layers.Dense(64, activation="relu", name="dense_1")(x)
        x = tf.keras.layers.Activation('tanh')(x)

        # x = tf.keras.layers.Dense(32, activation="relu", name="dense_2")(x)
        # x = tf.keras.layers.Dropout(0.6)(x)

        # x = tf.keras.layers.Dense(16, activation="relu", name="dense_3")(x)

        outputs = tf.keras.layers.Dense(5, name="predictions", bias_initializer=self.bias_init)(x)

        model = tf.keras.Model(inputs=inputs, outputs=outputs)

        model.compile(optimizer=tf.keras.optimizers.RMSprop(learning_rate=0.001),
                      loss=tf.keras.losses.SparseCategoricalCrossentropy(from_logits=True),
                      metrics=self.METRICS)
        return model

    def build_model(self, load_model: bool = False):
        """build model and run model training
        """
        if load_model:
            self.load_model(self.model_path)
        else:
            self.model = self.make_model_conv()
            self.model.evaluate(self.train_features, self.train_labels, batch_size=self.BATCH_SIZE, verbose=0)
            temp = tempfile.mkdtemp()
            initial_weights = os.path.join(temp, 'initial_weights/initial_weights')

            # print(temp);
            self.model.save_weights(initial_weights)
            #
            self.model = self.make_model_conv()
            self.model.load_weights(initial_weights)

#             weight_for_0 = (1 / self.train_class_counts['bi']) * (self.train_class_counts['total'] / 5.0)
#             weight_for_1 = (1 / self.train_class_counts['poo']) * (self.train_class_counts['total'] / 5.0)
#             weight_for_2 = (1 / self.train_class_counts['biased']) * (self.train_class_counts['total'] / 5.0)
#             weight_for_3 = (1 / self.train_class_counts['tb_arme']) * (self.train_class_counts['total'] / 5.0)
#             weight_for_4 = (1 / self.train_class_counts['zero_inflated']) * (self.train_class_counts['total'] / 5.0)
            
#             class_weights = {0: weight_for_0, 1: weight_for_1, 2: weight_for_2, 3: weight_for_3, 4: weight_for_4}

            # checkpoint_dir = os.path.dirname(self.checkpoint_path)
            # Create a callback that saves the model's weights
            cp_callback = tf.keras.callbacks.ModelCheckpoint(filepath=self.checkpoint_path, save_weights_only=True,
                                                             verbose=1)

            self.history = self.model.fit(
                self.train_features,
                self.train_labels,
                batch_size=self.BATCH_SIZE,
                epochs=self.EPOCHS,
                callbacks=[self.early_stopping, cp_callback],
                validation_data=(self.val_features, self.val_labels)
            )

#             print("bi count", self.train_class_counts['bi'])
#             print("poo count", self.train_class_counts['poo'])
#             print("biased count", self.train_class_counts['biased'])
#             print("tb_arme count", self.train_class_counts['tb_arme'])
#             print("zero_inflated count", self.train_class_counts['zero_inflated'])
#             print("total count", self.train_class_counts['total'])



        self.save_model(self.model_path)
        return

    def plot_metrics(self):
        """create line plots for training and validation loss
        """
        if self.history is not None:
            metrics = ['loss']
            for n, metric in enumerate(metrics):
                name = metric.replace("_", " ").capitalize()
                plt.subplot(2, 2, n + 1)
                plt.plot(self.history.epoch, self.history.history[metric], color=colors[0], label='Train')
                plt.plot(self.history.epoch, self.history.history['val_' + metric], color=colors[0], linestyle="--",
                         label='Val')
                plt.xlabel('Epoch')
                plt.ylabel(name)
                if metric == 'loss':
                    plt.ylim([0, plt.ylim()[1]])
                elif metric == 'auc':
                    plt.ylim([0.8, 1])
                else:
                    plt.ylim([0, 1])
                plt.legend()
        else:
            print("Training history not found. Possibly the model was not trained or a pretrained model was loaded " +
                  "to the system")

        return

    # def predict_classes(self, features: np.array, gene_names: list):
    #     """display model prediction for given set of input features and gene name list
    #
    #     Args:
    #         features (np.array): _description_
    #         gene_names (list): _description_
    #     """
    #     for index, feature in enumerate(features):
    #         result = self.model.predict(np.array([feature]))
    #         #             prediction = result
    #         prediction = tf.nn.softmax(result)
    #         print("This gene[{}] most likely belongs to class {} with a {:.2f} percent confidence.".format(
    #             gene_names[index], self.class_list[np.argmax(prediction)], 100 * np.max(prediction)))
    #     return

    def predict_classes(self, dataset):
        aggregate_df = dataset.fillna(0)

        if self.scaler is None:
            self.load_scaler(self.scaler_path)
        dataset_features = self.scaler.transform(np.array(aggregate_df))
        predict_classes = []
        for index, feature in enumerate(dataset_features):
            print(index, end=" ")
            result = self.model.predict(np.array([feature]), verbose=0)
            prediction = tf.nn.softmax(result)
            predict_classes.append(self.class_list[np.argmax(prediction)])
        aggregate_df.insert(0, 'predict_label', predict_classes)
        aggregate_df = aggregate_df[['predict_label']]
        aggregate_df.to_csv("data/all_prediction_only_2.csv", sep=";")

    def predict_anndata(self):
        aggregate_df = pd.DataFrame(self.adata.layers['aggregate'], columns=self.adata.var.index)
        aggregate_df = aggregate_df.transpose().fillna(0)

        if self.scaler is None:
            self.load_scaler(self.scaler_path)
        dataset_features = self.scaler.transform(np.array(aggregate_df))
        predict_classes = []
        for index, feature in enumerate(dataset_features):
            print(index, end=" ")
            result = self.model.predict(np.array([feature]), verbose=0)
            prediction = tf.nn.softmax(result)
            predict_classes.append(self.class_list[np.argmax(prediction)])

        self.adata.var['predicted_class'] = predict_classes

    def get_adata(self):
        return self.adata



    def evaluate_classes(self, prediction_set=None):
        """evaluate test dataset, and then save the predicted results in a column 'predict class'. Also displays
        the confusion matrix suing cross tab
        """
        if prediction_set is None:
            if self.test_dataset is None:
                print("No test dataset found. Use load_test_data() function to load the test dataset")
            predict_classes_test = []
            for index, feature in enumerate(self.test_features):
                result = self.model.predict(np.array([feature]), verbose=0)
                prediction = tf.nn.softmax(result)
#                 print(prediction)
                predict_classes_test.append(self.class_list[np.argmax(prediction)])

            self.test_dataset['predict_class'] = predict_classes_test

            print(pd.crosstab(self.test_dataset.class_names, self.test_dataset.predict_class))
        else:
            # for now prediction set is anndata
            allele_1 = pd.DataFrame(prediction_set.layers['spliced'], columns=prediction_set.var.index)
            allele_2 = pd.DataFrame(prediction_set.layers['unspliced'], columns=prediction_set.var.index)
            allele_1_T = allele_1.transpose(copy=True)
            allele_2_T = allele_2.transpose(copy=True)
            allele_T = (allele_1_T - allele_2_T) / (allele_1_T + allele_2_T)
            prediction_set.layers['aggregate'] = allele_T.transpose(copy=True)

            aggregate_layer = prediction_set.layers['aggregate']
            aggregate_layer_T = aggregate_layer.transpose()
            aggregate_df = pd.DataFrame(aggregate_layer_T, index=prediction_set.var.index)

            self.predict_classes(aggregate_df)

        return

    def summary(self):
        """shows network architecture summary
        """
        self.model.summary()
        return

    def evaluate(self):
        """evaluates test dataset
        """
        self.model.evaluate(x=self.test_features, y=self.test_labels)
        return
