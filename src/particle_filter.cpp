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

void ParticleFilter::init(double x, double y, double theta, double std[]) 
{
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
	default_random_engine gen;
	
	num_particles = 31;
	
	// This line creates a normal (Gaussian) distribution for x
	normal_distribution<double> dist_x(x, std[0]);
	// This line creates a normal (Gaussian) distribution for y
	normal_distribution<double> dist_y(y, std[1]);
	normal_distribution<double> dist_theta(theta, std[2]);
	
	for(int i = 0;i<num_particles;i++)
	{
	  Particle p;
	  p.id = i;
	  p.x = dist_x(gen);
	  p.y = dist_y(gen);
	  p.theta = dist_theta(gen);
	  p.weight = 1.0;
	  particles.push_back(p);
	  weights.push_back(p.weight);
	}
	is_initialized = true; 

}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) 
{
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
	default_random_engine gen;
	
	
	
   for (unsigned int i = 0; i < num_particles; i++)
   {
	  double p_x = particles[i].x;
	  double p_y = particles[i].y;
	  double p_theta = particles[i].theta;
	 
	  double pred_x;
	  double pred_y;
	  double pred_theta;
	  
	  //check for division by 0 for yaw_rate
	  if (fabs(yaw_rate) > 0.0001)
	  {
	  
	    pred_x = p_x + (velocity/yaw_rate) * (sin(p_theta + (yaw_rate * delta_t)) - sin(p_theta));
	    pred_y = p_y + (velocity/yaw_rate) * (cos(p_theta) - cos(p_theta + (yaw_rate * delta_t)));
	    pred_theta = p_theta + (yaw_rate * delta_t); 
	  }
	  else
	  {
	    pred_x = p_x + velocity * cos(p_theta) * delta_t;
	    pred_y = p_y + velocity * sin(p_theta) * delta_t;
	    pred_theta = p_theta;
	  }
	  
	  // This line creates a normal (Gaussian) distribution for x
	  normal_distribution<double> dist_x(pred_x, std_pos[0]);
	 // This line creates a normal (Gaussian) distribution for y
	  normal_distribution<double> dist_y(pred_y, std_pos[1]);
	  normal_distribution<double> dist_theta(pred_theta, std_pos[2]);
	  
	  particles[i].x = dist_x(gen);
	  particles[i].y = dist_y(gen);
	  particles[i].theta = dist_theta(gen);
	  
   }	  
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations)
{
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.
	// TODO complete

	//use nearest neighbor algorithm to associate observations in map co-ordinates to predicted landmarks  

	for (unsigned int i = 0; i < observations.size(); i++)
	{
	
	  //Set small_dist to infinity
	  double small_dist = std::numeric_limits<double>::infinity();
	  int nn_id = -1;
	  double obs_x = observations[i].x;
	  double obs_y = observations[i].y;

	  for (unsigned j = 0; j < predicted.size(); j++)
	  {
		double pred_x = predicted[j].x;
		double pred_y = predicted[j].y;
		int pred_id = predicted[j].id;
		double distance = dist(obs_x, obs_y, pred_x, pred_y);

		if (distance < small_dist)
		{
		  small_dist = distance;
		  nn_id = pred_id;
	    }
      }
	  observations[i].id = nn_id;
	}
}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		const std::vector<LandmarkObs> &observations, const Map &map_landmarks)
{
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
	
    for (unsigned int i = 0; i < num_particles; i++)
    {
       double p_x = particles[i].x;
       double p_y = particles[i].y;
       double p_theta = particles[i].theta;
  
       //Transform observations from vehicle co-ordinates to map co-ordinates.
       vector<LandmarkObs> p_observations;

       //Transform observations from vehicle's co-ordinates to map co-ordinates.
       for (unsigned int j = 0; j < observations.size(); j++)
       {
          LandmarkObs p_obs;
          p_obs.id = j;
          p_obs.x = p_x + (cos(p_theta) * observations[j].x) - (sin(p_theta) * observations[j].y);
          p_obs.y = p_y + (sin(p_theta) * observations[j].x) + (cos(p_theta) * observations[j].y);
          p_observations.push_back(p_obs);
       }

       //Filter map landmarks to keep ones within the sensor_range 
    
       vector<LandmarkObs> p_landmarks;
       for (unsigned int j = 0; j < map_landmarks.landmark_list.size(); j++)
       {
         double x_f = map_landmarks.landmark_list[j].x_f;
         double y_f = map_landmarks.landmark_list[j].y_f;
         int id_i = map_landmarks.landmark_list[j].id_i;
         if ((fabs((p_x - x_f)) <= sensor_range) && (fabs((p_y - y_f)) <= sensor_range)) {
           p_landmarks.push_back(LandmarkObs {id_i, x_f, y_f});
        }
     }

    //Associate particle observations to particle landmarks using nearest neighbor algorithm
    dataAssociation(p_landmarks, p_observations); //, sensor_range);

    //weight calculation of each particle using Multivariate Gaussian distribution.
   
    particles[i].weight = 1.0;

    double sd_x = std_landmark[0];
    double sd_y = std_landmark[1];
    double var_x = pow(sd_x, 2);
    double var_y = pow(sd_y, 2);
    double normalizer = (1.0/(2.0 * M_PI * sd_x * sd_y));
    
   
    for (unsigned int k = 0; k < p_observations.size(); k++)
    {
      double trans_obs_x = p_observations[k].x;
      double trans_obs_y = p_observations[k].y;
      double trans_obs_id = p_observations[k].id;
      double multi_prob = 1.0;

      for (unsigned l = 0; l < p_landmarks.size(); l++)
      {
        double pred_landmark_x = p_landmarks[l].x;
        double pred_landmark_y = p_landmarks[l].y;
        double pred_landmark_id = p_landmarks[l].id;

        if (trans_obs_id == pred_landmark_id)
        {
          multi_prob = normalizer * exp(-1.0 * ((pow((trans_obs_x - pred_landmark_x), 2)/(2.0 * var_x)) + (pow((trans_obs_y - pred_landmark_y), 2)/(2.0 * var_y))));
          particles[i].weight *= multi_prob;
        }
      }
    }
    
  }

   // ensure sum(weights) = 1
  double normal_factor = particles[0].weight;
  for (unsigned int i = 1; i < num_particles; i++)
  {
     normal_factor += particles[i].weight;
  }
      
  // Normalize the weights 
  for (int i = 0; i < particles.size(); i++)
  {
    particles[i].weight = particles[i].weight/normal_factor;
    weights[i] = particles[i].weight;
  }
}

void ParticleFilter::resample()
{
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
	
	vector<Particle> resampled_particles;

	// Create a generator to be used for generating random particle index and beta value
	default_random_engine gen;
	
	//Generate random particle index
	uniform_int_distribution<int> particle_index(0, num_particles - 1);
	
	int current_index = particle_index(gen);
	
	double beta = 0.0;
	
	double max_weight_2 = 2.0 * *max_element(weights.begin(), weights.end());
	
	for (int i = 0; i < particles.size(); i++)
	{
	  uniform_real_distribution<double> random_weight(0.0, max_weight_2);
	  beta += random_weight(gen);

	  while (beta > weights[current_index])
	  {
	    beta -= weights[current_index];
	    current_index = (current_index + 1) % num_particles;
	  }
	  resampled_particles.push_back(particles[current_index]);
	}
	particles = resampled_particles;

}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
                                     const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

  return particle;
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
