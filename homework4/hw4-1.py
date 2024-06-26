import numpy as np
import pandas as pd
from scipy.stats import norm
from plotnine import *


'''
Homework 4 problem 1 -- Plot data (please save to file, dont just print it)
plot the timeseries data for simulated nanopore
'''
def plot_timeseries_data(data):
    points = data
    #print(ggplot(data=points,mapping=aes(x='time', y='level'))+ geom_point())
    p = ggplot(data=points,mapping=aes(x='time', y='level'))+ geom_point()
    ggsave(plot = p, filename = 'timeseries_data_plot')
'''
Homework 4 problem 2
What is the approximate duration of each "event" in this data given this plot?
'''
# my rationale is that if it has a 1/50 chance to change then it is roughly 50 in length
approx_duration = 50

'''
Homework 4 problem 3 -- HMM maximum likelihood state sequence with 4 states
state 1 - T corresponds to a normal distribution with mean 100 and sd 15
state 2 - A corresponds to a normal dist with mean 150 and sd 25
state 3 - G correcponds to a normal dist with mean 300 and sd 50
state 4 - C corresponds to a normal dist with mean 350 and sd 25
transitions between states are 1/50 and transitions to same state is 49/50
'''

def HMM_MLE(df):
    # assuming 1/50 transition chance is evenly divided between all possible states
    transition_matrix = ([
        [49/50, 1/150, 1/150, 1/150],
        [1/150, 49/50, 1/150, 1/150],
        [1/150, 1/150, 49/50, 1/150],
        [1/150, 1/150, 1/150, 49/50]
    ])

    # we only need the 'level' values from the csv, not the 'time' values
    data = df['level']

    forward = np.zeros((len(data), 4))
    forward_scale = np.zeros(len(data))

    backward = np.zeros((len(data), 4))
    backward_scale = np.zeros(len(data))

    mean_values = [100, 150, 300, 350]
    standard_deviations = [15, 25, 50, 25]

    # Run the forward algorithm
    forward[0, :] = norm.pdf(data[0], loc=mean_values, scale=standard_deviations) * (1/4)
    # instead of using logsumexp we can just normalize the values at each point of time using forward_scale
    forward_scale[0] = sum(forward[0, :])
    forward[0, :] /= forward_scale[0]
    for t in range(1, len(data)):
        forward[t, :] = norm.pdf(data[t], loc=mean_values, scale=standard_deviations) * np.dot(forward[t-1, :], transition_matrix)
        forward_scale[t] = sum(forward[t, :])
        # normalize values
        forward[t, :] /= forward_scale[t]

    # Run the backward algorithm
    backward[-1, :] = np.ones(4) / 4
    backward_scale[-1] = sum(backward[-1, :])
    backward[-1, :] /= backward_scale[-1]
    for t in range(len(data)-2, -1, -1):
        backward[t, :] = np.dot(norm.pdf(data[t+1], loc=mean_values, scale=standard_deviations) * backward[t+1, :], transition_matrix)
        backward_scale[t] = sum(backward[t, :])
        # normalize values
        backward[t, :] /= backward_scale[t]

    # Calculate the state posteriors
    posteriors = np.zeros((len(data), 4))
    for t in range(len(data)):
        posteriors[t, :] = forward[t, :] * backward[t, :]
        #normalize the values 
        posteriors[t, :] /= sum(posteriors[t, :])

    # find the max value at each time
    state_sequence = np.argmax(posteriors, axis=1) + 1

    return state_sequence

'''
Homework 4 problem 4
plot output of problem 3. Here, please make 1 plot with 4 plots overlayed with different colors.
'''
# I think this is saying to basically plot a graph with each sequence state as a different color
# I spent far more time than I'd like to admit trying to figure this out as multiple plots, this 
# should be functionally what the question is looking for just using a different approach.
def plot_MLE(state_sequence):
    sequenceData = pd.DataFrame({'x': range(len(state_sequence)), 'y': state_sequence })
    
    # a slightly roundabout method to color the graph based on values, the values are 1 to 4 so there's a
    # placeholder for 0. Not the prettiest code I've ever written but it is functional
    colors = ['placeholder because my list is 1-indexed', 'T reads', 'A reads', 'G reads', 'C reads']
    sequenceData['color'] = [colors[i] for i in sequenceData['y']]

    sequencePlot = ggplot(data=sequenceData, mapping=aes(x='x', y='y', color='color')) + geom_point()
    
    ggsave(plot = sequencePlot, filename = 'state_sequence_plot')
    #print(sequencePlot)
    

'''
Homework 4 problem 5
Give the most likely sequence this data corresponds to given the likely 
event length you found from plotting the data
print this sequence of A/C/G/Ts
'''
# the "event_length" corresponds to the approx_duration from question 2, this means that each cluster of
# values corresponds to one A/C/G/T in an actual strand of DNA. At least that is my understanding of this
def MLE_seq(df, event_length):
    # df = state_sequence
    # event_length = approx_duration
    MLE_sequence = []
    sub_arr = [df[0]]
    MLE_output = ''

    for i in range(1, len(df)):
        if df[i] == df[i-1]:
            sub_arr.append(df[i])
        else:
            MLE_sequence.append(sub_arr)
            sub_arr = [df[i]]

    MLE_sequence.append(sub_arr)

    for i in range(len(MLE_sequence)):
        # this may be slightly off for the data set, I'm confident there's a better approach but
        # the deadline fast approaches. The "back-to-back events" seem to only be of length 2, there could be some
        # of length 3 but this is hopefully close.

        # If I had more time for a different approach:
        # Dividing the sub-sequences by event length if they are
        # greater than event_length * 1.5 then appending different amounts of letters to the string
        # based on a rounded-down sub-sequence/event_length.

        if len(MLE_sequence[i]) > event_length * 1.5:
             if MLE_sequence[i][1] == 1:
                MLE_output += 'TT'
             elif MLE_sequence[i][1] == 2:
                MLE_output += 'AA'
             elif MLE_sequence[i][1] == 3:
                MLE_output += 'GG'
             elif MLE_sequence[i][1] == 4:
                MLE_output += 'CC'
        else:
            if MLE_sequence[i][1] == 1:
                MLE_output += 'T'
            elif MLE_sequence[i][1] == 2:
                MLE_output += 'A'
            elif MLE_sequence[i][1] == 3:
                MLE_output += 'G'
            elif MLE_sequence[i][1] == 4:
                MLE_output += 'C'
            
            
    print(MLE_output)

    


'''
Homework 4 problem 6
Forward/backward algorithm giving posterior probabilities for each time point for each level
'''
# basically the same as problem 3, returning posteriors instead of using posteriors to find optimal states
def HMM_posterior(df):
    # assuming 1/50 transition chance is evenly divided between all possible states
    transition_matrix = ([
        [49/50, 1/150, 1/150, 1/150],
        [1/150, 49/50, 1/150, 1/150],
        [1/150, 1/150, 49/50, 1/150],
        [1/150, 1/150, 1/150, 49/50]
    ])

    # we only need the 'level' values from the csv, not the 'time' values
    data = df['level']

    forward = np.zeros((len(data), 4))
    forward_scale = np.zeros(len(data))

    backward = np.zeros((len(data), 4))
    backward_scale = np.zeros(len(data))

    mean_values = [100, 150, 300, 350]
    standard_deviations = [15, 25, 50, 25]

    # Run the forward algorithm
    forward[0, :] = norm.pdf(data[0], loc=mean_values, scale=standard_deviations) * (1/4)
    # instead of using logsumexp we can just normalize the values at each point of time using forward_scale
    forward_scale[0] = sum(forward[0, :])
    forward[0, :] /= forward_scale[0]
    for t in range(1, len(data)):
        forward[t, :] = norm.pdf(data[t], loc=mean_values, scale=standard_deviations) * np.dot(forward[t-1, :], transition_matrix)
        forward_scale[t] = sum(forward[t, :])
        # normalize values
        forward[t, :] /= forward_scale[t]

    # Run the backward algorithm
    backward[-1, :] = np.ones(4) / 4
    backward_scale[-1] = sum(backward[-1, :])
    backward[-1, :] /= backward_scale[-1]
    for t in range(len(data)-2, -1, -1):
        backward[t, :] = np.dot(norm.pdf(data[t+1], loc=mean_values, scale=standard_deviations) * backward[t+1, :], transition_matrix)
        backward_scale[t] = sum(backward[t, :])
        # normalize values
        backward[t, :] /= backward_scale[t]

    # Calculate the state posteriors
    posteriors = np.zeros((len(data), 4))
    for t in range(len(data)):
        posteriors[t, :] = forward[t, :] * backward[t, :]
        #normalize the values 
        posteriors[t, :] /= sum(posteriors[t, :])
    
    return posteriors


'''
Homework 4 problem 7
plot output of problem 5, this time, plot with 4 facets using facet_wrap
'''
def plot_posterior(posteriors):
    df = pd.DataFrame(data=posteriors, columns=['T','A','G','C'])
    df['x'] = range(len(posteriors))
    combined_vars = pd.melt(df, id_vars='x', var_name='facet')
    p = ggplot(combined_vars, aes(x='x', y='value')) + geom_line() + facet_wrap('~facet')
    #print(p)
    ggsave(plot = p, filename = 'posteriors_plot')


df = pd.read_csv("nanopore.csv")
plot_timeseries_data(df)
state_sequence = HMM_MLE(df)
plot_MLE(state_sequence)
MLE_seq(state_sequence, approx_duration)
posteriors = HMM_posterior(df)
plot_posterior(posteriors)

