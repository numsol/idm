#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/* select precision. */
#if (defined SINGLE_PRECISION)
  typedef float decimal;
#elif (defined DOUBLE_PRECISION)
  typedef double decimal;
#elif (defined QUADRUPLE_PRECISION)
  typedef __float128 decimal;
#else
  typedef float decimal;
#endif /* PRECISION */

/* units of length are meters.
   units of time are seconds. */

/* car properties. */
typedef struct car_t
{
  decimal position_x;
  decimal velocity_x;
  decimal desired_velocity; 
  decimal desired_time_headway;
  decimal max_acceleration;
  decimal desired_deceleration;
  decimal acceleration_exponent;
  decimal jam_distance;
  decimal vehicle_length;
  decimal net_distance;
  decimal approaching_rate;
} car;

decimal
position_derivative (car o[])
{
  return (o->velocity_x);
}

decimal
velocity_derivative (car o[])
{
  decimal term_a = exp (o->acceleration_exponent * \
    log (o->velocity_x / o->desired_velocity));

  decimal term_b_numerator_a = o->jam_distance + \
    (o->velocity_x * o->desired_time_headway);

  decimal term_b_numerator_b = \
    (o->velocity_x * o->approaching_rate) / \
    (2 * sqrt (o->max_acceleration * o->desired_deceleration));

  decimal term_b_numerator_sum = term_b_numerator_a + \
    term_b_numerator_b;

  decimal term_b_dividend = term_b_numerator_sum / \
    o->net_distance;

  decimal term_b = exp (2 * log (term_b_dividend));

  return (o->max_acceleration * (1 - term_a - term_b));
}

int
main (int argc, char **argv)
{
  int n_cars = 20;
  decimal n_steps = 50000;
  decimal current_time = 0.0;
  decimal time_step = 0.1;
  decimal road_length = 400.0;

  /* loop over the cars to set the defaults in an array
     of car parameter structures. */
  int i;
  car car_object[n_cars];
  for (i = 0; i < n_cars; i++)
  {
    /* distribute the cars over the length of the road
       equally. we don't know the net distace or the 
       approaching rate yet, so we leave these alone
       for now. */
    car_object[i].position_x = \
      ((decimal)i / (decimal)n_cars) * road_length;
    car_object[i].velocity_x = 20.0;
    car_object[i].desired_velocity = 30.0; 
    car_object[i].desired_time_headway = 1.5;
    car_object[i].max_acceleration = 1.0;
    car_object[i].desired_deceleration = 3.0;
    car_object[i].acceleration_exponent = 4.0;
    car_object[i].jam_distance = 2.0;
    car_object[i].vehicle_length = 5.0;
  }

  /* positions and velocities files. */
  FILE *positions, *velocities;
  positions = fopen ("positions.dat", "a");
  velocities = fopen ("velocities.dat", "a");

  /* the phase point is set up so we can now start the loop. */
  int j, next_car;
  decimal k1, k2, k3, k4;
  decimal m1, m2, m3, m4;
  decimal velocity_orig;
  for (i = 0; i < n_steps + 1; i++)
  {
    /* print the timestep to file. */
    fprintf (positions, "%20.6f", (float)current_time);
    fprintf (velocities, "%20.6f", (float)current_time);

    /* build the current net distance and approaching rate. */
    for (j = 0; j < n_cars; j++)
    {
      /* print the positions and velocities to file. */
      fprintf (positions, "%20.6f", (float)car_object[j].position_x);
      fprintf (velocities, "%20.6f", (float)car_object[j].velocity_x);

      /* define the next car with periodic boundary conditions
         so the first car in the loop is seen by the last. */
      next_car = (j + 1) % n_cars;

      car_object[j].net_distance = \
        car_object[next_car].position_x - \
        car_object[j].position_x - \
        car_object[next_car].vehicle_length;

      /* if the distance is less than 0, then the next car has 
         been moved across the periodic boundary. */
      if (car_object[j].net_distance < 0)
        car_object[j].net_distance += road_length;

      car_object[j].approaching_rate = \
        car_object[j].velocity_x - \
        car_object[next_car].velocity_x;
    }

      /* begin runge-kutta order four. */
    for (j = 0; j < n_cars; j++)
    {
      /* store the original values of the velocity
         so it can be modified from and restored 
         during the rk procedure. the positions do
         do not appear in the derivatives and 
         therefore need not be treated here. */
      velocity_orig = car_object[j].velocity_x;

      k1 = position_derivative (&car_object[j]);
      m1 = velocity_derivative (&car_object[j]);

      car_object[j].velocity_x = \
        velocity_orig + (0.5 * m1 * time_step);

      k2 = position_derivative (&car_object[j]);
      m2 = velocity_derivative (&car_object[j]);

      car_object[j].velocity_x = \
        velocity_orig + (0.5 * m2 * time_step);

      k3 = position_derivative (&car_object[j]);
      m3 = velocity_derivative (&car_object[j]);

      car_object[j].velocity_x = \
        velocity_orig + (m3 * time_step);

      k4 = position_derivative (&car_object[j]);
      m4 = velocity_derivative (&car_object[j]);

      /* finally update the variables. */
      current_time += time_step;
      car_object[j].position_x += \
        ((k1 + (2 * k2) + (2 * k3) + k4)) * (time_step / 6);
      car_object[j].velocity_x = velocity_orig + \
        ((m1 + (2 * m2) + (2 * m3) + m4)) * (time_step / 6);

      /* if the car has moved beyond the end of the road, 
         move it back to the beginning. */
      if (car_object[j].position_x > road_length)
        car_object[j].position_x -= road_length;
    }

    /* we're done with this step's output so print
       a new line in output files. */
    fprintf (positions, "\n");
    fprintf (velocities, "\n");
  }

  /* close file pointers. */
  fclose (positions);
  fclose (velocities);

  return EXIT_SUCCESS;
}
