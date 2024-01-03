from typing import Literal, Tuple
import numpy as np
from scipy.optimize import minimize_scalar


class AbstractSpringlock:
    # Helper variables to convert units
    _unit_conversions = {
        'length': {'m': 1.000, 'cm': 1.000e-2, 'mm': 1.000e-3, 'ft': 0.3048, 'in': 0.0254, 'banana': 0.1842},
        'mass': {'kg': 1.000, 'g': 1.000e-3, 'lb': 0.4536, 'oz': 0.02835},
        'force': {'N': 1.000, 'lbf': 4.448, 'ozf': 0.2780},
        'spring constant': {'N/m': 1.000, 'lbf/in': 175.1, 'N/mm': 1.000e3, 'lbf/ft': 14.59},
        'speed': {'m/s': 1.000, 'ft/s': 0.3048, 'mph': 0.4470},
        'torque': {'N*m': 1.000, 'lbf*ft': 1.356, 'ozf*ft': 0.08474, 'lbf*in': 0.1130, 'ozf*in': 0.007062},
        'angle': {'rad': 1.000, 'deg': 0.01745},
        'angular velocity': {'rad/s': 1.000, 'deg/s': 0.01745, 'rpm': 0.1047},
        'energy': {'J': 1.000, 'N*m': 1.000, 'lbf*ft': 1.356, 'ozf*ft': 0.08474, 'lbf*in': 0.1130, 'ozf*in': 0.007062},
        'pressure': {'N/m^2': 1.000, 'Pa': 1.000, 'lbf/in^2': 6.895e3, 'psi': 6.895e3, 'bar': 1.000e5, 'atm': 1.013e5, 'MPa': 1.000e6},
        'area': {'m^2': 1.000, 'cm^2': 1.000e-4, 'mm^2': 1.000e-6, 'ft^2': 0.09290, 'in^2': 6.452e-4},
        'volume': {'m^3': 1.000, 'cm^3': 1.000e-6, 'mm^3': 1.000e-9, 'ft^3': 2.832e-2, 'in^3': 1.638e-5, 'l': 1.000e-3, 'ml': 1.000e-6}
    }

    _unit_types = {
        'length': Literal['m', 'cm', 'mm', 'ft', 'in', 'banana'],
        'mass': Literal['kg', 'g', 'lb', 'oz'],
        'force': Literal['N', 'lbf', 'ozf'],
        'spring constant': Literal['N/m', 'lbf/in', 'N/mm', 'lbf/ft'],
        'speed': Literal['m/s', 'ft/s', 'mph'],
        'torque': Literal['N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in'],
        'angle': Literal['rad', 'deg'],
        'angular velocity': Literal['rad/s', 'deg/s', 'rpm'],
        'energy': Literal['J', 'N*m', 'lbf*ft', 'ozf*ft', 'lbf*in', 'ozf*in'],
        'pressure': Literal['N/m^2', 'Pa', 'lbf/in^2', 'psi', 'bar', 'atm', 'MPa'],
        'area': Literal['m^2', 'cm^2', 'mm^2', 'ft^2', 'in^2'],
        'volume': Literal['m^3', 'cm^3', 'mm^3', 'ft^3', 'in^3', 'l', 'ml']
    }

    def __init__(self, name: str,
                    weapon_mass: np.float64, weapon_radius: np.float64, weapon_arm_length: np.float64, weapon_arm_mass: np.float64, spring_arm_length: np.float64, spring_arm_mass: np.float64,
                 *, length_units: _unit_types['length'] = 'm', mass_units: _unit_types['mass'] = 'kg', force_units: _unit_types['force'] = 'N',
                    spring_constant_units: _unit_types['spring constant'] = 'N/m', speed_units: _unit_types['speed'] = 'm/s', torque_units: _unit_types['torque'] = 'N*m',
                    angle_units: _unit_types['angle'] = 'rad', angular_velocity_units: _unit_types['angular velocity'] = 'rad/s', energy_units: _unit_types['energy'] = 'J',
                    pressure_units:_unit_types['pressure'] = 'N/m^2', area_units: _unit_types['area'] = 'm^2', volume_units: _unit_types['volume'] = 'm^3'):
        
        self.name = name

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
        self.pressure_units = pressure_units
        self.area_units = area_units
        self.volume_units = volume_units

        self.weapon_mass = weapon_mass
        self.weapon_radius = weapon_radius
        self.weapon_arm_length = weapon_arm_length
        self.weapon_arm_mass = weapon_arm_mass
        self.spring_arm_length = spring_arm_length
        self.spring_arm_mass = spring_arm_mass

        # Class variables that will be calculated and memoized as needed
        self._moment_of_inertia = None
        self._max_torque_theta = None
        self._max_radial_force_theta = None
        self._max_tangential_force_theta = None
        self._max_spring_speed_theta = None


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


    ###################################
    ###### Unit conversion stuff ######
    ###################################
        

    @property
    def length_units(self) -> _unit_types['length']:
        return self._units['length']
    
    @length_units.setter
    def length_units(self, value: _unit_types['length']):
        self._units['length'] = value


    @property
    def mass_units(self) -> _unit_types['mass']:
        return self._units['mass']
    
    @mass_units.setter
    def mass_units(self, value: _unit_types['mass']):
        self._units['mass'] = value


    @property
    def force_units(self) -> _unit_types['force']:
        return self._units['force']
    
    @force_units.setter
    def force_units(self, value: _unit_types['force']):
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


    @property
    def pressure_units(self) -> Literal['N/m^2', 'Pa', 'lbf/in^2', 'psi', 'bar', 'atm', 'MPa']:
        return self._units['pressure']
    
    @pressure_units.setter
    def pressure_units(self, value: Literal['N/m^2', 'Pa', 'lbf/in^2', 'psi', 'bar', 'atm', 'MPa']):
        self._units['pressure'] = value


    @property
    def area_units(self) -> Literal['m^2', 'cm^2', 'mm^2', 'ft^2', 'in^2']:
        return self._units['area']
    
    @area_units.setter
    def area_units(self, value: Literal['m^2', 'cm^2', 'mm^2', 'ft^2', 'in^2']):
        self._units['area'] = value

    
    @property
    def volume_units(self) -> Literal['m^3', 'cm^3', 'mm^3', 'ft^3', 'in^3', 'l', 'ml']:
        return self._units['volume']
    
    @volume_units.setter
    def volume_units(self, value: Literal['m^3', 'cm^3', 'mm^3', 'ft^3', 'in^3', 'l', 'ml']):
        self._units['volume'] = value



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
    def max_tangential_force_theta(self) -> np.float64:
        if self._max_tangential_force_theta is None:
            self._max_tangential_force_theta = self._calculate_max_tangential_force_theta()
        return self._convert_units(self._max_tangential_force_theta, 'angle', False)
    
    @property
    def max_spring_speed_theta(self) -> np.float64:
        if self._max_spring_speed_theta is None:
            self._max_spring_speed_theta = self._calculate_max_spring_speed_theta()
        return self._convert_units(self._max_spring_speed_theta, 'angle', False)
    
    @property
    def max_torque(self) -> np.float64:
        return self.torque(self.max_torque_theta)
    
    @property
    def max_radial_force(self) -> np.float64:
        return self.radial_force(self.max_radial_force_theta)
    
    @property
    def max_tangential_force(self) -> np.float64:
        return self.tangential_force(self.max_tangential_force_theta)
    
    @property
    def max_spring_speed(self) -> np.float64:
        return self.spring_speed(self.max_spring_speed_theta)
    

    @property
    def max_tip_speed(self) -> np.float64:
        return self.tip_speed(0)   

    @property
    def max_potential_energy(self) -> np.float64:
        return self._convert_units(self._max_potential_energy, 'energy', False)
    
    @property
    def _max_potential_energy(self) -> np.float64:
        return self._potential_energy(np.pi) - self._potential_energy(0)
    
    

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
        raise NotImplementedError()
    

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
        raise NotImplementedError()
    

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
        raise NotImplementedError()


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
        energy_delta = np.abs(max_energy - current_energy)

        return np.sqrt(2 * energy_delta / self._moment_of_inertia)
    

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
        raise NotImplementedError()
    

    def tangential_force(self, arm_angle: np.float64) -> np.float64:
        """Calculates the tangential force on the arm at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The tangential force on the arm.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_tangential_force = self._tangential_force(base_angle)
        return self._convert_units(base_tangential_force, 'force', False)
    
    def _tangential_force(self, arm_angle: np.float64) -> np.float64:
        raise NotImplementedError()
    

    def spring_speed(self, arm_angle: np.float64) -> np.float64:
        """Calculates the linear velocity of the spring/spring piston at a given angle theta
            when the springlock is released from the fully loaded position
        
        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The linear velocity of the spring/spring piston.
        """
        
        base_angle = self._convert_units(arm_angle, 'angle')
        base_speed = self._spring_speed(base_angle)
        return self._convert_units(base_speed, 'speed', False)
    
    def _spring_speed(self, arm_angle: np.float64) -> np.float64:
        raise NotImplementedError()
    


    ############################
    ###### Helper Methods ######
    ############################


    def _convert_units(self, value: np.float64, unit_type: Literal['length', 'mass', 'force', 'spring constant', 'speed', 'torque', 'angle', 'angular velocity'], to_base: bool = True) -> np.float64:
        units_scale = AbstractSpringlock._unit_conversions[unit_type][self._units[unit_type]]
        return value * units_scale if to_base else value / units_scale
    

    def _clear_state(self):
        self._moment_of_inertia = None
        self._max_torque_theta = None
        self._max_radial_force_theta = None
        self._max_tangential_force_theta = None


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
        

    def _calculate_max_tangential_force_theta(self):
        """Calculates the maximum tangential force on the arm.

        Returns:
            np.float64: The maximum tangential force on the arm.
        """            
        optimize_results = minimize_scalar(lambda x: -abs(self._tangential_force(x)), bounds=(0, np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum tangential force of springlock system')
        

    def _calculate_max_spring_speed_theta(self):
        """Calculates the maximum linear velocity of the spring/spring piston.

        Returns:
            np.float64: The maximum linear velocity of the spring/spring piston.
        """            
        optimize_results = minimize_scalar(lambda x: -abs(self._spring_speed(x)), bounds=(0, np.pi), method='bounded')

        if optimize_results.success:
            return optimize_results.x
        else:
            raise RuntimeError('Unable to find maximum spring speed of springlock system')



##############################################################################
######                       GAS SPRING SPRINGLOCK                      ######
##############################################################################



class GasSpringlock(AbstractSpringlock):
    def __init__(self, name: str,
                    weapon_mass: np.float64, weapon_radius: np.float64, weapon_arm_length: np.float64, weapon_arm_mass: np.float64, spring_arm_length: np.float64, spring_arm_mass: np.float64,
                    spring_initial_length: np.float64, spring_final_length: np.float64, spring_initial_force: np.float64, spring_final_force: np.float64, spring_initial_pressure: np.float64,
                *, length_units: AbstractSpringlock._unit_types['length'] = 'm', mass_units: AbstractSpringlock._unit_types['mass'] = 'kg', force_units: AbstractSpringlock._unit_types['force'] = 'N',
                    spring_constant_units: AbstractSpringlock._unit_types['spring constant'] = 'N/m', speed_units: AbstractSpringlock._unit_types['speed'] = 'm/s', torque_units: AbstractSpringlock._unit_types['torque'] = 'N*m',
                    angle_units: AbstractSpringlock._unit_types['angle'] = 'rad', angular_velocity_units: AbstractSpringlock._unit_types['angular velocity'] = 'rad/s', energy_units: AbstractSpringlock._unit_types['energy'] = 'J',
                    pressure_units: AbstractSpringlock._unit_types['pressure'] = 'N/m^2', area_units: AbstractSpringlock._unit_types['area'] = 'm^2', volume_units: AbstractSpringlock._unit_types['volume'] = 'm^3'):
        
        # Initialize the unit stuff in the base class
        super().__init__(
            name=name,
            weapon_mass=weapon_mass,
            weapon_radius=weapon_radius,
            weapon_arm_length=weapon_arm_length,
            weapon_arm_mass=weapon_arm_mass,
            spring_arm_length=spring_arm_length,
            spring_arm_mass=spring_arm_mass,
            length_units=length_units,
            mass_units=mass_units,
            force_units=force_units,
            spring_constant_units=spring_constant_units,
            speed_units=speed_units,
            torque_units=torque_units,
            angle_units=angle_units,
            angular_velocity_units=angular_velocity_units,
            energy_units=energy_units,
            pressure_units=pressure_units,
            area_units=area_units,
            volume_units=volume_units
            )
        
        self.spring_initial_length = spring_initial_length
        self.spring_final_length = spring_final_length
        self.spring_initial_force = spring_initial_force
        self.spring_final_force = spring_final_force
        self.spring_initial_pressure = spring_initial_pressure


    ###################################
    ###### Getter/Setter Methods ######
    ###################################
        
    
    @property
    def spring_initial_length(self) -> np.float64:
        return self._convert_units(self._spring_initial_length, 'length', False)
    
    @spring_initial_length.setter
    def spring_initial_length(self, value: np.float64):
        self._spring_initial_length = self._convert_units(value, 'length')
        self._clear_state()

    
    @property
    def spring_final_length(self) -> np.float64:
        return self._convert_units(self._spring_final_length, 'length', False)
    
    @spring_final_length.setter
    def spring_final_length(self, value: np.float64):
        self._spring_final_length = self._convert_units(value, 'length')
        self._clear_state()


    @property
    def spring_initial_force(self) -> np.float64:
        return self._convert_units(self._spring_initial_force, 'force', False)
    
    @spring_initial_force.setter
    def spring_initial_force(self, value: np.float64):
        self._spring_initial_force = self._convert_units(value, 'force')
        self._clear_state()


    @property
    def spring_final_force(self) -> np.float64:
        return self._convert_units(self._spring_final_force, 'force', False)
    
    @spring_final_force.setter
    def spring_final_force(self, value: np.float64):
        self._spring_final_force = self._convert_units(value, 'force')
        self._clear_state()


    @property
    def spring_initial_pressure(self) -> np.float64:
        return self._convert_units(self._spring_initial_pressure, 'pressure', False)
    
    @spring_initial_pressure.setter
    def spring_initial_pressure(self, value: np.float64):
        self._spring_initial_pressure = self._convert_units(value, 'pressure')
        self._clear_state()



    ###################################
    ###### Calculated Properties ######
    ###################################
        

    @property
    def piston_area(self) -> np.float64:
        return self._convert_units(self._piston_area, 'area', False)
    
    @property
    def _piston_area(self) -> np.float64:
        return self._spring_initial_force / self._spring_initial_pressure
    

    @property
    def spring_final_pressure(self) -> np.float64:
        return self._convert_units(self._spring_final_pressure, 'pressure', False)
    
    @property
    def _spring_final_pressure(self) -> np.float64:
        return self._spring_final_force / self._piston_area
    

    @property
    def spring_volume_change(self) -> np.float64:
        return self._convert_units(self._spring_volume_change, 'volume', False)
    
    @property
    def _spring_volume_change(self) -> np.float64:
        return (self._spring_final_length - self._spring_initial_length) * self._piston_area
    

    @property
    def spring_initial_volume(self) -> np.float64:
        return self._convert_units(self._spring_initial_volume, 'volume', False)
    
    @property
    def _spring_initial_volume(self) -> np.float64:
        return (self._spring_final_pressure * self._spring_volume_change) / (self._spring_initial_pressure - self._spring_final_pressure)


    @property
    def spring_final_volume(self) -> np.float64:
        return self._convert_units(self._spring_final_volume, 'volume', False)
    
    @property
    def _spring_final_volume(self) -> np.float64:
        return self._spring_initial_volume + self._spring_volume_change

    #############################
    ###### Physics Methods ######
    #############################


    def _spring_length(self, arm_angle: np.float64) -> np.float64:
        length_a = self._spring_initial_length - self._spring_arm_length
        length_b = self._spring_arm_length
        return np.sqrt(length_a**2 + length_b**2 + 2 * length_a * length_b * np.cos(arm_angle))
    
    
    def spring_volume(self, arm_angle: np.float64) -> np.float64:
        """Calculates the volume of the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The volume of the spring.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_volume = self._spring_volume(base_angle)
        return self._convert_units(base_volume, 'volume', False)
    
    def _spring_volume(self, arm_angle: np.float64) -> np.float64:
        return self._spring_initial_volume + self._piston_area * (self._spring_length(arm_angle) - self._spring_initial_length)
    

    def spring_pressure(self, arm_angle: np.float64) -> np.float64:
        """Calculates the pressure of the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The pressure of the spring.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_pressure = self._spring_pressure(base_angle)
        return self._convert_units(base_pressure, 'pressure', False)

    def _spring_pressure(self, arm_angle: np.float64) -> np.float64:
        return self._spring_initial_pressure * self._spring_initial_volume / self._spring_volume(arm_angle)
    

    def spring_force(self, arm_angle: np.float64) -> np.float64:
        """Calculates the force of the spring at a given arm angle.

        Args:
            arm_angle (np.float64): The angle of the arm

        Returns:
            np.float64: The force of the spring.
        """
        base_angle = self._convert_units(arm_angle, 'angle')
        base_force = self._spring_force(base_angle)
        return self._convert_units(base_force, 'force', False)
    
    def _spring_force(self, arm_angle: np.float64) -> np.float64:
        return self._spring_pressure(arm_angle) * self._piston_area
    

    def _potential_energy(self, arm_angle: np.float64) -> np.float64:
        log_term = self._spring_initial_volume + self._piston_area * (self._spring_length(arm_angle) - self._spring_initial_length)
        log_term /= self._spring_initial_volume
        log_term  = np.abs(log_term)
        return self._spring_initial_pressure * self._spring_initial_volume * np.log(log_term)
    

    def _torque(self, arm_angle: np.float64) -> np.float64:
        leg_a = self._spring_initial_length - self._spring_arm_length
        leg_b = self._spring_arm_length
    
        numerator = self._spring_initial_pressure * self._spring_initial_volume * self._piston_area
        numerator *= leg_a * leg_b * np.sin(arm_angle)

        denominator = self._spring_initial_volume + self._piston_area * (self._spring_length(arm_angle) - self._spring_initial_length)
        denominator *= self._spring_length(arm_angle)

        return numerator / denominator


    def _radial_force(self, arm_angle: np.float64) -> np.float64:
        return np.sqrt(self._spring_force(arm_angle)**2 - self._tangential_force(arm_angle)**2)
    

    def _tangential_force(self, arm_angle: np.float64) -> np.float64:
        return self._torque(arm_angle) / self._weapon_arm_length
    

    def _spring_speed(self, arm_angle: np.float64) -> np.float64:
        arm_a = self._spring_initial_length - self._spring_arm_length
        arm_b = self._spring_arm_length
        angular_velocity = self._angular_velocity(arm_angle)

        return arm_a * arm_b * angular_velocity * np.sin(arm_angle) / self._spring_length(arm_angle)
    

