import numpy as np

def twobody(t: float, x: np.ndarray, mu: float = 3.986004418E5) -> np.ndarray:
    """Returns the two body dynamics of an object in orbit around Earth.

    Args:
        t (float): epoch in seconds
        x (np.ndarray): state = [r, v] in units [km, km/s]
        mu (float, optional): Gravitational Parameter. Defaults to 3.986004418E5.

    Returns:
        np.ndarray: rate of change of state
    """
    
    r = x[:3]
    v = x[3:]
    
    return np.concatenate([v, -mu * r / np.linalg.norm(r)**3])

def twobody_STM(t: float, x: np.ndarray, mu: float = 3.986004418E5) -> np.ndarray:
    """Returns the two body dynamics of an object in orbit around Earth. Also propagates the State Transition Matrix

    Args:
        t (float): epoch in seconds
        x (np.ndarray): state = [r, v, Phi_flattened] in units [km, km/s]
        mu (float, optional): Gravitational Parameter. Defaults to 3.986004418E5

    Returns:
        np.ndarray: rates of change in states and STM
    """
    r = x[:3]
    v = x[3:6]
    
    Phi = x[6:].reshape((6,6))
    
    Z3 = np.zeros((3,3))
    I3 = np.eye(3)
    rnorm = np.linalg.norm(r)
    
    A = np.block(
        [[Z3, I3], 
         [3*mu/rnorm**5 * (r @ r.T) - mu/rnorm**3 * I3, Z3]]
        )
    
    Phidot = A @ Phi
    
    return np.concatenate([v, -mu * r / rnorm**3, Phidot.flatten()])

def inertial_range(r: np.ndarray, R: float, lat: float, LST: float) -> np.ndarray:
    """Calculates slant range vector resolved in inertial frame

    Args:
        r (np.ndarray): _description_
        R (float): _description_
        lat (float): _description_
        LST (float): _description_

    Returns:
        np.ndarray: slant range vector
    """
    rotation = np.array([np.cos(lat)*np.cos(LST), np.cos(lat)*np.sin(LST), np.sin(lat)])
    
    rho = r - R * rotation
    
    return rho

def observer_range(rho: np.ndarray, lat: float, LST: float):
    """_summary_

    Args:
        rho (np.ndarray): _description_
        lat (float): _description_
        LST (float): _description_
    """
    
    rho_obsv = np.zeros_like(rho)
    
    R2 = np.array([
        [np.cos(lat), 0, np.sin(lat)],
        [0, 1, 0],
        [-np.siN(lat), 0, np.cos(lat)]
    ])
    
    for k in range(len(LST)):
        LST_k = LST[k]
        R1 = np.array([
            [np.cos(LST_k), np.sin(LST_k), 0],
            [-np.sin(LST_k), np.cos(LST_k), 0],
            [0, 0, 1]
        ])

def topocentric():
    pass

def add_noise():
    pass
