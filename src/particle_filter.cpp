/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).

	default_random_engine gen;

	num_particles = 100;

	// Creates a normal (Gaussian) distribution for x,y and theta.
	normal_distribution<double> dist_x(x, std[0]);
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);

    for (int i = 1; i <= num_particles; i++){
        Particle p = {};    
        // assign values to p
        p.id = i;
        p.x = dist_x(gen);
        p.y = dist_y(gen);
        p.theta = dist_theta(gen);
        p.weight = 1.0;

        // add p to your particles vector
        particles.push_back(p);
    }

    is_initialized = true;

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/

	default_random_engine gen;

	  // define normal distributions for sensor noise with zero mean and standard deviationa as std_pos
    normal_distribution<double> N_x_noise(0, std_pos[0]);
    normal_distribution<double> N_y_noise(0, std_pos[1]);
    normal_distribution<double> N_theta_noise(0, std_pos[2]);

    for (int i = 0; i < num_particles; i++) {

        if(yaw_rate < 0.0001)
        {
            particles[i].x += velocity * delta_t * cos(particles[i].theta);
            particles[i].y += velocity * delta_t * sin(particles[i].theta);
            /* Theta remains same since yaw rate is almost zero */
        }
        else
        {
            particles[i].x += velocity / yaw_rate * (sin(particles[i].theta + yaw_rate*delta_t) - sin(particles[i].theta));
            particles[i].y += velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate*delta_t));
            particles[i].theta += yaw_rate * delta_t;

        }
        
        /* Add noise */
        particles[i].x += N_x_noise(gen);
        particles[i].y += N_y_noise(gen);
        particles[i].theta += N_theta_noise(gen);
    }


}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

    /* Loop though all the observed landmarks */
    for (unsigned int i = 0; i < observations.size(); i++) {
    
       LandmarkObs obs_pos = observations[i];

       /* Initialize minimum distance to max value */
       double min_dist = numeric_limits<double>::max();

       /* Initialize ID of the nearest observed landmark */
       int id = -1;
    
       /* For each landmark observation, find the nearest predicted observation */
       for (unsigned int j = 0; j < predicted.size(); j++) {
      
           LandmarkObs pred_pos = predicted[j];
      
           /* Compute the distance between current/predicted landmarks */
           double cur_dist = dist(obs_pos.x, obs_pos.y, pred_pos.x, pred_pos.y);

           /* Identify the landmark nearest to the predicted landmark observation */
           if (cur_dist < min_dist) {
              min_dist = cur_dist;
              id = pred_pos.id;
            }
        }
    /* set the observation's id to the nearest predicted landmark's id */
    observations[i].id = id;
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33
	//   http://planning.cs.uiuc.edu/node99.html

    vector<LandmarkObs> landmarks_in_range;
    vector<LandmarkObs> map_observations;

    // Loop through each particle
    for (int i = 0; i < num_particles; ++i){

        vector<LandmarkObs> transformed_obs;
        // For each observation
        for (int j = 0; j < observations.size(); ++j){
            // Transform the observation point (from vehicle coordinates to map coordinates)
            double trans_obs_x, trans_obs_y;
            trans_obs_x = observations[j].x * cos(particles[i].theta) - observations[j].y * sin(particles[i].theta) + particles[i].x;
            trans_obs_y = observations[j].x * sin(particles[i].theta) + observations[j].y * cos(particles[i].theta) + particles[i].y;

            transformed_obs.push_back(LandmarkObs{ observations[j].id, trans_obs_x, trans_obs_x});
        }

        /* Find map landmarks that are within sensor range */
        for (int j = 0;  j < map_landmarks.landmark_list.size(); j++){
            int landmark_id = map_landmarks.landmark_list[j].id_i;
            float landmark_x = map_landmarks.landmark_list[j].x_f;
            float landmark_y = map_landmarks.landmark_list[j].y_f;

            float dist_x = landmark_x - particles[i].x;
            float dist_y = landmark_y - particles[i].y;

            float range = sqrt(dist_x + dist_y);
            if(range < sensor_range){

                 landmarks_in_range.push_back(LandmarkObs{landmark_id, landmark_x, landmark_y});
            }
        }

        /* Perform data association */
        dataAssociation(landmarks_in_range, transformed_obs);

    }      
    
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
