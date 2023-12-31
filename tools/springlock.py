from typing import Literal, Tuple
import numpy as np
from scipy.optimize import minimize_scalar


class Springlock:
    # Helper variables to convert units
    _unit_conversions = {
        'length': {'m': 1.000, 'cm': 1.000e-2, 'mm': 1.000e-3, 'ft': 0.3048, 'in': 0.0254},
        'mass': {'kg': 1.000, 'g': 1.000e-3, 'lb': 0.4536, 'oz': 0.02835},
        'force': {'N': 1.000, 'lbf': 4.448, 'ozf': 0.2780},
        'spring constant': {'N/m': 1.000, 'lbf/in': 175.1, 'N/mm': 1.000e3, 'lbf/ft': 14.59},
        'speed': {'m/s': 1.000, 'ft/s': 0.3048, 'mph': 0.4470},
        'torque': {'N*m': 1.000, 'lbf*ft': 1.356, 'ozf*ft': 0.08474, 'lbf*in': 0.1130, 'ozf*in': 0.007062},
        'angle': {'rad': 1.000, 'deg': 0.01745},
        'angular velocity': {'rad/s': 1.000, 'deg/s': 0.01745, 'rpm': 0.1047},
        'energy': {'J': 1.000}
    }


    def __init__(self, weapon_mass: np.float64, weapon_radius: np.float64, weapon_arm_length: np.float64, weapon_arm_mass: np.float64, spring_arm_length: np.float64, 
                    spring_arm_mass: np.float64, spring_constant: np.float64, spring_resting_length: np.float64, spring_stroke_length: np.float64,
                 *, length_units: Literal['m', 'cm', 'mm', 'ft', 'in'] = 'm', mass_units: Literal['kg', 'g', 'lb', 'oz'] = 'kg', force_units: Literal['N', 'lbf', 'ozf'] = 'N',
                    spring_constant_units: Literal['N/m', 'lbf/in', 'N/mm', 'lbf/ft'] = 'N/m', speed_units: Literal['m/s', 'ft/s', 'mph'] = 'm/s',
                    torque_units: Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in'] = 'N*m', angle_units: Literal['rad', 'deg'] = 'rad',
                    angular_velocity_units: Literal['rad/s', 'deg/s', 'rpm'] = 'rad/s', energy_units: Literal['J'] = 'J'):
        
        """Creates a new Springlock object with the specified parameters.

        Args:
            weapon_radius (np.float64): The distance from the base of the arm to the center of the weapon.
            weapon_mass (np.float64): The mass of the weapon.
            spring_radius (np.float64): The distance from the base of the arm to the mounting point of the spring.
            spring_constant (np.float64): The spring constant of the springs. If multiple springs are used in parallel, this should be the sum of the sprting constants.
            spring_resting_angle (np.float64): The angle at which the spring when it is not under any tension or compression.
            spring_mounting_point (Tuple[np.float64, np.float64]): The x and y coordinates of the fixed mounting point of the spring, relative to the base of the arm
        """
        self._units = {}
        self.length_units = length_units
        self.mass_units = mass_units
        self.force_units = force_units
        self.spring_constant_units = spring_constant_units
        self.speed_units = speed_units
        self.torque_units = torque_units
        self.angle_units = angle_units
        self.angular_velocity_units = angular_velocity_units
        self.energy_units = energy_units

        self.weapon_mass = weapon_mass
        self.weapon_radius = weapon_radius
        self.weapon_arm_length = weapon_arm_length
        self.weapon_arm_mass = weapon_arm_mass
        self.spring_arm_length = spring_arm_length
        self.spring_arm_mass = spring_arm_mass
        self.spring_constant = spring_constant
        self.spring_resting_length = spring_resting_length
        self.spring_stroke_length = spring_stroke_length

        # Class variables that will be calculated and memoized as needed
        self._moment_of_inertia = None
        self._max_torque_theta = None
        self._max_radial_force_theta = None
        self._max_angular_force_theta = None



    ###################################
    ###### Getter/Setter Methods ######
    ###################################
    
    
    @property
    def weapon_mass(self) -> np.float64:
        return self._convert_units(self._weapon_mass, 'mass', False)
    
    @weapon_mass.setter
    def weapon_mass(self, value: np.float64):
        self._weapon_mass = self._convert_units(value, 'mass')
        self._clear_state()


    @property
    def weapon_radius(self) -> np.float64:
        return self._convert_units(self._weapon_radius, 'length', False)
    
    @weapon_radius.setter
    def weapon_radius(self, value: np.float64):
        self._weapon_radius = self._convert_units(value, 'length')
        self._clear_state()


    @property
    def weapon_arm_length(self) -> np.float64:
        return self._convert_units(self._weapon_arm_length, 'length', False)

    @weapon_arm_length.setter
    def weapon_arm_length(self, value: np.float64):
        self._weapon_arm_length = self._convert_units(value, 'length')
        self._clear_state()


    @property
    def weapon_arm_mass(self) -> np.float64:
        return self._convert_units(self._weapon_arm_mass, 'mass', False)
    
    @weapon_arm_mass.setter
    def weapon_arm_mass(self, value: np.float64):
        self._weapon_arm_mass = self._convert_units(value, 'mass')
        self._clear_state()


    @property
    def spring_arm_length(self) -> np.float64:
        return self._convert_units(self._spring_arm_length, 'length', False)
    
    @spring_arm_length.setter
    def spring_arm_length(self, value: np.float64):
        self._spring_arm_length = self._convert_units(value, 'length')
        self._clear_state()


    @property
    def spring_arm_mass(self) -> np.float64:
        return self._convert_units(self._spring_arm_mass, 'mass', False)
    
    @spring_arm_mass.setter
    def spring_arm_mass(self, value: np.float64):
        self._spring_arm_mass = self._convert_units(value, 'mass')
        self._clear_state()

    
    @property
    def spring_constant(self) -> np.float64:
        return self._convert_units(self._spring_constant, 'spring constant', False)
    
    @spring_constant.setter
    def spring_constant(self, value: np.float64):
        self._spring_constant = self._convert_units(value, 'spring constant')
        self._clear_state()


    @property
    def spring_resting_length(self) -> np.float64:
        return self._convert_units(self._spring_resting_length, 'length', False)
    
    @spring_resting_length.setter
    def spring_resting_length(self, value: np.float64):
        self._spring_resting_length = self._convert_units(value, 'length')
        self._clear_state()


    @property
    def spring_stroke_length(self) -> np.float64:
        return self._convert_units(self._spring_stroke_length, 'length', False)
    
    @spring_stroke_length.setter
    def spring_stroke_length(self, value: np.float64):     
        if self._convert_units(value, 'length') < self._spring_arm_length * 2:
            raise ValueError('Spring stroke length must be greater than twice the spring arm length')
        
        self._spring_stroke_length = self._convert_units(value, 'length')
        self._clear_state()



    ###################################
    ###### Calculated Properties ######
    ###################################


    @property
    def moment_of_inertia(self) -> np.float64:
        if self._moment_of_inertia is None:
            self._moment_of_inertia = self._calculate_moment_of_inertia()
        mass_scale = self._convert_units(1, 'mass', False)
        length_scale = self._convert_units(1, 'length', False)
        return self._moment_of_inertia / mass_scale * length_scale**2
    
    @property
    def max_torque_theta(self) -> np.float64:
        if self._max_torque_theta is None:
            self._max_torque_theta = self._calculate_max_torque_theta()
        return self._convert_units(self._max_torque_theta, 'angle', False)
    
    @property
    def max_radial_force_theta(self) -> np.float64:
        if self._max_radial_force_theta is None:
            self._max_radial_force_theta = self._calculate_max_radial_force_theta()
        return self._convert_units(self._max_radial_force_theta, 'angle', False)
    
    @property
    def max_angular_force_theta(self) -> np.float64:
        if self._max_angular_force_theta is None:
            self._max_angular_force_theta = self._calculate_max_angular_force_theta()
        return self._convert_units(self._max_angular_force_theta, 'angle', False)
    
    @property
    def max_torque(self) -> np.float64:
        return self.torque(self.max_torque_theta)
    
    @property
    def max_radial_force(self) -> np.float64:
        return self.radial_force(self.max_radial_force_theta)
    
    @property
    def max_angular_force(self) -> np.float64:
        return self.angular_force(self.max_angular_force_theta)
    
    @property
    def max_tip_speed(self) -> np.float64:
        return self.tip_speed(0)
    

    
    ###################################
    ###### Unit conversion stuff ######
    ###################################
        

    @property
    def length_units(self) -> Literal['m', 'cm', 'mm', 'ft', 'in']:
        return self._units['length']
    
    @length_units.setter
    def length_units(self, value: Literal['m', 'cm', 'mm', 'ft', 'in']):
        self._units['length'] = value


    @property
    def mass_units(self) -> Literal['kg', 'g', 'lb', 'oz']:
        return self._units['mass']
    
    @mass_units.setter
    def mass_units(self, value: Literal['kg', 'g', 'lb', 'oz']):
        self._units['mass'] = value


    @property
    def force_units(self) -> Literal['N', 'lbf', 'ozf']:
        return self._units['force']
    
    @force_units.setter
    def force_units(self, value: Literal['N', 'lbf', 'ozf']):
        self._units['force'] = value

    @property
    def spring_constant_units(self) -> Literal['N/m', 'lbf/in', 'N/mm', 'lbf/ft']:
        return self._units['spring constant']
    
    @spring_constant_units.setter
    def spring_constant_units(self, value: Literal['N/m', 'lbf/in', 'N/mm', 'lbf/ft']):
        self._units['spring constant'] = value


    @property
    def speed_units(self) -> Literal['m/s', 'ft/s', 'mph']:
        return self._units['speed']
    
    @speed_units.setter
    def speed_units(self, value: Literal['m/s', 'ft/s', 'mph']):
        self._units['speed'] = value


    @property
    def torque_units(self) -> Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in']:
        return self._units['torque']
    
    @torque_units.setter
    def torque_units(self, value: Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in']):
        self._units['torque'] = value


    @property
    def angle_units(self) -> Literal['rad', 'deg']:
        return self._units['angle']
    
    @angle_units.setter
    def angle_units(self, value: Literal['rad', 'deg']):
        self._units['angle'] = value


    @property
    def angular_velocity_units(self) -> Literal['rad/s', 'deg/s', 'rpm']:
        return self._units['angular velocity']
    
    @angular_velocity_units.setter
    def angular_velocity_units(self, value: Literal['rad/s', 'deg/s', 'rpm']):
        self._units['angular velocity'] = value


    @property
    def energy_units(self) -> Literal['J']:
        return 'J'
    
    @energy_units.setter
    def energy_units(self, value: Literal['J']):
        self._units['energy'] = value
    


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
        base_angle = self._convert_units(arm_angle, 'angle')
        base_length = self._spring_length(base_angle)
        return self._convert_units(base_length, 'length', False)

    def _spring_length(self, arm_angle: np.float64) -> np.float64:
        length_a = self._spring_arm_length
        length_b = self._spring_arm_length + self._spring_minimum_length
        return np.sqrt(length_a**2 + length_b**2 - 2 * length_a * length_b * np.cos(arm_angle))
    

    def potential_energy(self, arm_angle: np.float64) -> np.float64:
        """Calculates the energy stored in the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The energy stored in the spring.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_energy = self._potential_energy(base_angle)
        return self._convert_units(base_energy, 'energy', False)
    
    def _potential_energy(self, arm_angle: np.float64) -> np.float64:
        return 0.5 * self._spring_constant * (self._spring_length(arm_angle) - self._spring_minimum_length)**2
    

    def torque(self, arm_angle: np.float64) -> np.float64:
        """Calculates the torque on the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The torque on the arm.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_torque = self._torque(base_angle)
        return self._convert_units(base_torque, 'torque', False)
    
    def _torque(self, arm_angle: np.float64) -> np.float64:
        leg_a = self._spring_arm_length
        leg_b = self._spring_arm_length + self._spring_minimum_length
        spring_extension = self._spring_length(arm_angle) - self._spring_resting_length
        return -self._spring_constant * leg_a * leg_b * np.sin(arm_angle) * spring_extension / self._spring_length(arm_angle)
    

    def angular_velocity(self, arm_angle: np.float64) -> np.float64:
        """Calculates the angular velocity of the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The angular velocity of the arm.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_angular_velocity = self._angular_velocity(base_angle)
        return self._convert_units(base_angular_velocity, 'angular velocity', False)
    
    def _angular_velocity(self, arm_angle: np.float64) -> np.float64:
        if self._moment_of_inertia is None:
            self._moment_of_inertia = self._calculate_moment_of_inertia()

        max_energy = self._potential_energy(np.pi)
        current_energy = self._potential_energy(arm_angle)

        return np.sqrt(2 * (max_energy - current_energy) / self._moment_of_inertia)
    
    
    def tip_speed(self, arm_angle: np.float64) -> np.float64:
        """Calculates the speed of the tip of the arm at a given arm angle when released from the maximum potential energy position.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The velocity of the arm.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_speed = self._tip_speed(base_angle)
        return self._convert_units(base_speed, 'speed', False)
    
    def _tip_speed(self, arm_angle: np.float64) -> np.float64:
        return self._angular_velocity(arm_angle) * self._weapon_arm_length
    

    def radial_force(self, arm_angle: np.float64) -> np.float64:
        """Calculates the radial force on the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The radial force on the arm.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_radial_force = self._radial_force(base_angle)
        return self._convert_units(base_radial_force, 'force', False)
    
    def _radial_force(self, arm_angle: np.float64) -> np.float64:
        leg_a = self._spring_arm_length
        leg_b = self._spring_length(arm_angle)
        leg_c = self._spring_arm_length + self._spring_minimum_length
        spring_extension = self._spring_length(arm_angle) - self._spring_resting_length
        return self._spring_constant * spring_extension * (leg_c**2 - leg_a**2 - leg_b**2) / (2 * leg_a * leg_b)
    

    def angular_force(self, arm_angle: np.float64) -> np.float64:
        """Calculates the angular force on the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The angular force on the arm.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_angular_force = self._angular_force(base_angle)
        return self._convert_units(base_angular_force, 'force', False)
    
    def _angular_force(self, arm_angle: np.float64) -> np.float64:
        leg = self._spring_arm_length + self._spring_minimum_length
        spring_extension = self._spring_length(arm_angle) - self._spring_resting_length
        return -self._spring_constant * leg * np.sin(arm_angle) * spring_extension / self._spring_length(arm_angle)



    ############################
    ###### Helper Methods ######
    ############################


    def _convert_units(self, value: np.float64, unit_type: Literal['length', 'mass', 'force', 'spring constant', 'speed', 'torque', 'angle', 'angular velocity'], to_base: bool = True) -> np.float64:
        units_scale = Springlock._unit_conversions[unit_type][self._units[unit_type]]
        return value * units_scale if to_base else value / units_scale


    def _clear_state(self):
        self._moment_of_inertia = None
        self._max_torque_theta = None
        self._max_radial_force_theta = None
        self._max_angular_force_theta = None

    def _calculate_moment_of_inertia(self) -> np.float64:
        """Calculates the moment of inertia of the arm.

        Returns:
            np.float64: The moment of inertia of the arm.
        """
        weapon_arm_moi = 1/3 * self._weapon_arm_mass * self._weapon_arm_length**2
        spring_arm_moi = 1/3 * self._spring_arm_mass * self._spring_arm_length**2
        weapon_moi = 1/2 * self._weapon_mass * self._weapon_radius**2 + self._weapon_mass * self._weapon_arm_length**2

        return weapon_arm_moi + spring_arm_moi + weapon_moi


    def _calculate_max_torque_theta(self):
        """Calculates the maximum torque on the arm.

        Returns:
            np.float64: The maximum torque on the arm.
        """            
        optimize_results = minimize_scalar(lambda x: -abs(self._torque(x)), bounds=(0, np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum torque of springlock system')
        

    def _calculate_max_radial_force_theta(self):
        """Calculates the maximum radial force on the arm.

        Returns:
            np.float64: The maximum radial force on the arm.
        """            
        optimize_results = minimize_scalar(lambda x: -abs(self._radial_force(x)), bounds=(0, np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum radial force of springlock system')
        

    def _calculate_max_angular_force_theta(self):
        """Calculates the maximum angular force on the arm.

        Returns:
            np.float64: The maximum angular force on the arm.
        """            
        optimize_results = minimize_scalar(lambda x: -abs(self._angular_force(x)), bounds=(0, np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum angular force of springlock system')

