from typing import Literal, Tuple
import numpy as np
from scipy.optimize import minimize_scalar


class Springlock:
    # Helper variables to convert units
    _length_unit_conversions = {
        'm': 1.000,
        'cm': 1.000e-2,
        'mm': 1.000e-3,
        'ft': 0.3048,
        'in': 0.0254
    }

    _mass_unit_conversions = {
        'kg': 1.000,
        'g': 1.000e-3,
        'lb': 0.4536,
        'oz': 0.02835
    }

    _spring_constant_unit_conversions = {
        'N/m': 1.000,
        'lbf/in': 175.1
    }

    _speed_unit_conversions = {
        'm/s': 1.000,
        'ft/s': 0.3048,
        'mph': 0.4470
    }

    _torque_unit_conversions = {
        'N*m': 1.000,
        'lbf*ft': 1.356,
        'ozf*ft': 0.08474,
        'lbf*in': 0.1130,
        'ozf*in': 0.007062
    }

    _angle_unit_conversions = {
        'rad': 1.000,
        'deg': 0.01745
    }


    def __init__(self, weapon_radius: np.float64, weapon_mass: np.float64, spring_radius: np.float64, spring_constant: np.float64, spring_resting_angle: np.float64, spring_mounting_point: Tuple[np.float64, np.float64],
                 *, length_units: Literal['m', 'cm', 'mm', 'ft', 'in'] = 'm', mass_units: Literal['kg', 'g', 'lb', 'oz'] = 'kg', spring_constant_units: Literal['N/m', 'lbf/in'] = 'N/m',
                    speed_units: Literal['m/s', 'ft/s', 'mph'] = 'm/s', torque_units: Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in'] = 'N*m', angle_units: Literal['rad', 'deg'] = 'rad'):
        
        """Creates a new Springlock object with the specified parameters.

        Args:
            weapon_radius (np.float64): The distance from the base of the arm to the center of the weapon.
            weapon_mass (np.float64): The mass of the weapon.
            spring_radius (np.float64): The distance from the base of the arm to the mounting point of the spring.
            spring_constant (np.float64): The spring constant of the springs. If multiple springs are used in parallel, this should be the sum of the sprting constants.
            spring_resting_angle (np.float64): The angle at which the spring when it is not under any tension or compression.
            spring_mounting_point (Tuple[np.float64, np.float64]): The x and y coordinates of the fixed mounting point of the spring, relative to the base of the arm
        """
        

        self._length_scale = Springlock._length_unit_conversions[length_units]
        self._mass_scale = Springlock._mass_unit_conversions[mass_units]
        self._spring_constant_scale = Springlock._spring_constant_unit_conversions[spring_constant_units]
        self._speed_scale = Springlock._speed_unit_conversions[speed_units]
        self._torque_scale = Springlock._torque_unit_conversions[torque_units]
        self._angle_scale = Springlock._angle_unit_conversions[angle_units]

        # Class variables that will not be changed
        self._weapon_radius = weapon_radius * self._length_scale
        self._weapon_mass = weapon_mass * self._mass_scale
        self._spring_radius = spring_radius * self._length_scale
        self._spring_constant = spring_constant * self._spring_constant_scale
        self._spring_resting_angle = spring_resting_angle * self._angle_scale
        self._spring_mounting_point = (spring_mounting_point[0] * self._length_scale, spring_mounting_point[1] * self._length_scale)

        self._length_units = length_units
        self._mass_units = mass_units
        self._spring_constant_units = spring_constant_units
        self._speed_units = speed_units
        self._torque_units = torque_units
        self._angle_units = angle_units

        # Class variables that will be calculated and memoized as needed
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
        self._max_spring_theta = None
        self._min_spring_theta = None


    ###################################
    ###### Getter/Setter Methods ######
    ###################################
    
    @property
    def weapon_radius(self) -> np.float64:
        return self._weapon_radius / self._length_scale
    
    @weapon_radius.setter
    def weapon_radius(self, weapon_radius: np.float64):
        self._weapon_radius = weapon_radius * self._length_scale
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def weapon_mass(self) -> np.float64:
        return self._weapon_mass / self._mass_scale
    
    @weapon_mass.setter
    def weapon_mass(self, weapon_mass: np.float64):
        self._weapon_mass = weapon_mass * self._mass_scale
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def spring_radius(self) -> np.float64:
        return self._spring_radius / self._length_scale
    
    @spring_radius.setter
    def spring_radius(self, spring_radius: np.float64):
        self._spring_radius = spring_radius * self._length_scale
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def spring_constant(self) -> np.float64:
        return self._spring_constant / self._spring_constant_scale
    
    @spring_constant.setter
    def spring_constant(self, spring_constant: np.float64):
        self._spring_constant = spring_constant * self._spring_constant_scale
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def spring_resting_angle(self) -> np.float64:
        return self._spring_resting_angle / self._angle_scale
    
    @spring_resting_angle.setter
    def spring_resting_angle(self, spring_resting_angle: np.float64):
        self._spring_resting_angle = spring_resting_angle * self._angle_scale
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def spring_mounting_point(self) -> Tuple[np.float64, np.float64]:
        return (self._spring_mounting_point[0] / self._length_scale, self._spring_mounting_point[1] / self._length_scale)
    
    @spring_mounting_point.setter
    def spring_mounting_point(self, spring_mounting_point: Tuple[np.float64, np.float64]):
        self._spring_mounting_point = (spring_mounting_point[0] * self._length_scale, spring_mounting_point[1] * self._length_scale)
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def spring_mounting_point_x(self) -> np.float64:
        return self._spring_mounting_point[0] / self._length_scale
    
    @spring_mounting_point_x.setter
    def spring_mounting_point_x(self, spring_mounting_point_x: np.float64):
        self._spring_mounting_point = (spring_mounting_point_x * self._length_scale, self._spring_mounting_point[1])
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def spring_mounting_point_y(self) -> np.float64:
        return self._spring_mounting_point[1] / self._length_scale
    
    @spring_mounting_point_y.setter
    def spring_mounting_point_y(self, spring_mounting_point_y: np.float64):
        self._spring_mounting_point = (self._spring_mounting_point[0], spring_mounting_point_y * self._length_scale)
        self._max_theta = None
        self._min_theta = None
        self._max_torque_theta = None
    
    @property
    def tip_to_spring(self) -> np.float64:
        return self.weapon_radius - self.spring_radius
    
    @property
    def max_theta(self) -> np.float64:
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()
        return self._max_theta / self._angle_scale
    
    @property
    def release_angle(self) -> np.float64:
        return self.max_theta
    
    @property
    def max_potential_energy(self) -> np.float64:
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()
        return self._total_potential_energy(self._max_theta)
    
    @property
    def min_theta(self) -> np.float64:
        if self._min_theta is None:
            self._min_theta = self._calculate_min_theta()
        return self._min_theta / self._angle_scale
    
    @property
    def impact_angle(self) -> np.float64:
        return self.min_theta
    
    @property
    def min_potential_energy(self) -> np.float64:
        if self._min_theta is None:
            self._min_theta = self._calculate_min_theta()
        return self._total_potential_energy(self._min_theta)
    
    @property
    def max_speed(self) -> np.float64:
        return (np.sqrt(2 * (self.max_potential_energy - self.min_potential_energy) / self._weapon_mass)) / self._speed_scale
    
    @property
    def max_torque_theta(self) -> np.float64:
        if self._max_torque_theta is None:
            self._max_torque_theta = self._calculate_max_torque_theta()
        return self._max_torque_theta / self._angle_scale
    
    @property
    def max_torque(self) -> np.float64:
        if self._max_torque_theta is None:
            self._max_torque_theta = self._calculate_max_torque_theta()
        return np.abs(self._torque(self._max_torque_theta)) / self._torque_scale
    
    @property
    def resting_spring_length(self) -> np.float64:
        return self._spring_length(self._spring_resting_angle) / self._length_scale
    
    @property
    def max_spring_length(self) -> np.float64:
        if self._max_spring_theta is None:
            self._max_spring_theta = self._calculate_max_spring_theta()

        return self._spring_length(self._max_spring_theta) / self._length_scale
    
    @property
    def min_spring_length(self) -> np.float64:
        if self._min_spring_theta is None:
            self._min_spring_theta = self._calculate_min_spring_theta()

        return self._spring_length(self._min_spring_theta) / self._length_scale

    

    ###################################
    ###### Unit conversion stuff ######
    ###################################

    @property
    def length_units(self) -> Literal['m', 'cm', 'mm', 'ft', 'in']:
        return self._length_units
    
    @length_units.setter
    def length_units(self, value: Literal['m', 'cm', 'mm', 'ft', 'in']):
        self._length_units = value
        self._length_scale = Springlock._length_unit_conversions[value]

    @property
    def mass_units(self) -> Literal['kg', 'g', 'lb', 'oz']:
        return self._mass_units
    
    @mass_units.setter
    def mass_units(self, value: Literal['kg', 'g', 'lb', 'oz']):
        self._mass_units = value
        self._mass_scale = Springlock._mass_unit_conversions[value]

    @property
    def spring_constant_units(self) -> Literal['N/m', 'lbf/in']:
        return self._spring_constant_units
    
    @spring_constant_units.setter
    def spring_constant_units(self, value: Literal['N/m', 'lbf/in']):
        self._spring_constant_units = value
        self._spring_constant_scale = Springlock._spring_constant_unit_conversions[value]


    @property
    def speed_units(self) -> Literal['m/s', 'ft/s', 'mph']:
        return self._speed_units
    
    @speed_units.setter
    def speed_units(self, value: Literal['m/s', 'ft/s', 'mph']):
        self._speed_units = value
        self._speed_scale = Springlock._speed_unit_conversions[value]

    @property
    def torque_units(self) -> Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in']:
        return self._torque_units
    
    @torque_units.setter
    def torque_units(self, value: Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in']):
        self._torque_units = value
        self._torque_scale = Springlock._torque_unit_conversions[value]

    @property
    def angle_units(self) -> Literal['rad', 'deg']:
        return self._angle_units
    
    @angle_units.setter
    def angle_units(self, value: Literal['rad', 'deg']):
        self._angle_units = value
        self._angle_scale = Springlock._angle_unit_conversions[value]

    @property
    def energy_units(self) -> Literal['J']:
        return 'J'
    

    #############################
    ###### Physics Methods ######
    #############################

    def spring_length(self, arm_angle: np.float64) -> np.float64:
        """Calculates the length of the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The length of the spring.
        """
        return self._spring_length(arm_angle * self._angle_scale) / self._length_scale

    def _spring_length(self, arm_angle: np.float64) -> np.float64:
        """Calculates the length of the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The length of the spring.
        """
        return np.sqrt((self._spring_radius * np.cos(arm_angle) - self._spring_mounting_point[0])**2 + (self._spring_radius * np.sin(arm_angle) - self._spring_mounting_point[1])**2)
    

    def spring_potential_energy(self, arm_angle: np.float64) -> np.float64:
        """Calculates the energy stored in the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The energy stored in the spring.
        """
        return self._spring_potential_energy(arm_angle * self._angle_scale)
    
    def _spring_potential_energy(self, arm_angle: np.float64) -> np.float64:
        return 0.5 * self._spring_constant * (self._spring_length(arm_angle) - self._spring_length(self._spring_resting_angle))**2
    

    def gravitational_potential_energy(self, arm_angle: np.float64) -> np.float64:
        """Calculates the gravitational potential energy of the weapon at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The gravitational potential energy of the weapon.
        """
        return self._gravitational_potential_energy(arm_angle * self._angle_scale)
    
    def _gravitational_potential_energy(self, arm_angle: np.float64) -> np.float64:
        return self._weapon_mass * 9.81 * self._weapon_radius * (np.sin(arm_angle) - np.sin(self._spring_resting_angle))
    

    def total_potential_energy(self, arm_angle: np.float64) -> np.float64:
        """Calculates the total potential energy of the system at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The total potential energy of the system.
        """
        return self._total_potential_energy(arm_angle * self._angle_scale)
    
    def _total_potential_energy(self, arm_angle: np.float64) -> np.float64:
        return self._spring_potential_energy(arm_angle) + self._gravitational_potential_energy(arm_angle)
    

    def torque(self, arm_angle: np.float64) -> np.float64:
        """Calculates the torque on the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The torque on the arm.
        """
        return self._torque(arm_angle * self._angle_scale) / self._torque_scale
    
    def _torque(self, arm_angle: np.float64) -> np.float64:
        spring_torque = self._spring_constant * self._spring_radius * (self._spring_mounting_point[0] * np.sin(arm_angle) - self._spring_mounting_point[1] * np.cos(arm_angle)) * (self._spring_length(arm_angle) - self._spring_length(self._spring_resting_angle)) / self._spring_length(arm_angle)
        gravity_torque = self._weapon_mass * 9.81 * self._weapon_radius * np.cos(arm_angle)
        return -(spring_torque + gravity_torque)
    

    def angular_velocity(self, arm_angle: np.float64) -> np.float64:
        """Calculates the angular velocity of the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The angular velocity of the arm.
        """
        return self._angular_velocity(arm_angle * self._angle_scale) / self._angle_scale
    
    def _angular_velocity(self, arm_angle: np.float64) -> np.float64:
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()

        max_energy = self._total_potential_energy(self._max_theta)
        current_energy = self._total_potential_energy(arm_angle)

        return np.sqrt(2 * (max_energy - current_energy) / (self._weapon_mass * self._weapon_radius**2))
    
    
    def speed(self, arm_angle: np.float64) -> np.float64:
        """Calculates the speed of the tip of the arm at a given arm angle when released from the maximum potential energy position.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The velocity of the arm.
        """
        return self._speed(arm_angle * self._angle_scale) / self._speed_scale
    
    def _speed(self, arm_angle: np.float64) -> np.float64:
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()

        max_energy = self._total_potential_energy(self._max_theta)
        current_energy = self._total_potential_energy(arm_angle)

        return np.sqrt(2 * (max_energy - current_energy) / self._weapon_mass)


    ############################
    ###### Helper Methods ######
    ############################


    def _calculate_max_theta(self) -> np.float64:
        """Calculates the angle of the arm at which the total potential energy of the system is maximized.

        Returns:
            np.float64: The maximum angle of the arm in radians.
        """
        optimize_results = minimize_scalar(lambda x: -self._total_potential_energy(x), bounds=(-np.pi, 2 * np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum potential energy of springlock system')
    

    def _calculate_min_theta(self) -> np.float64:
        """Calculates the angle of the arm at which the total potential energy of the system is minimized.

        Returns:
            np.float64: The minimum angle of the arm in radians.
        """
        optimize_results = minimize_scalar(lambda x: self._total_potential_energy(x), bounds=(-2 * np.pi, np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find minimum potential energy of springlock system')


    def _calculate_max_torque_theta(self) -> np.float64:
        """Calculates the maximum torque on the arm.

        Returns:
            np.float64: The maximum torque on the arm.
        """

        if self._min_theta is None:
            self._min_theta = self._calculate_min_theta()
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()

        lower_bound = min(self._min_theta, self._max_theta)
        upper_bound = max(self._min_theta, self._max_theta)
            
        optimize_results = minimize_scalar(lambda x: self._torque(x), bounds=(lower_bound, upper_bound), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum torque of springlock system')


    def _calculate_max_spring_theta(self) -> np.float64:
        """Calculates the angle of the arm at which the spring is at maximum extension.

        Returns:
            np.float64: The angle of the arm at which the spring is at maximum extension.
        """

        if self._min_theta is None:
            self._min_theta = self._calculate_min_theta()
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()

        lower_bound = min(self._min_theta, self._max_theta)
        upper_bound = max(self._min_theta, self._max_theta)

        optimize_results = minimize_scalar(lambda x: self._spring_length(x), bounds=(lower_bound, upper_bound), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum extension of springlock system')
        
        
    def _calculate_min_spring_theta(self) -> np.float64:
        """Calculates the angle of the arm at which the spring is at maximum compression.

        Returns:
            np.float64: The angle of the arm at which the spring is at maximum compression.
        """

        if self._min_theta is None:
            self._min_theta = self._calculate_min_theta()
        if self._max_theta is None:
            self._max_theta = self._calculate_max_theta()

        lower_bound = min(self._min_theta, self._max_theta)
        upper_bound = max(self._min_theta, self._max_theta)

        optimize_results = minimize_scalar(lambda x: -self._spring_length(x), bounds=(lower_bound, upper_bound), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find minimum extension of springlock system')