import numpy as np
from scipy import interpolate

import importlib

import constants as const
import grid

const = importlib.reload(const)
grid = importlib.reload(grid)


# %% interpolation
def lin_lin_interp(xp, yp, kind='linear', axis=-1):
    interp = interpolate.interp1d(xp, yp, kind=kind, axis=axis)
    def func(x): return interp(x)
    return func


def lin_log_interp(xp, yp, kind='linear', axis=-1):
    log_yp = np.log10(yp)
    interp = interpolate.interp1d(xp, log_yp, kind=kind, axis=axis)
    def func(x): return np.power(10.0, interp(x))
    return func


def log_lin_interp(xp, yp, kind='linear', axis=-1):
    log_xp = np.log10(xp)
    interp = interpolate.interp1d(log_xp, yp, kind=kind, axis=axis)
    def func(x): return interp(np.log10(x))
    return func


def log_log_interp(xp, yp, kind='linear', axis=-1):
    log_xp = np.log10(xp)
    log_yp = np.log10(yp)
    interp = interpolate.interp1d(log_xp, log_yp, kind=kind, axis=axis)
    def func(x): return np.power(10.0, interp(np.log10(x)))
    return func


def lin_log_log_interp_2d(xp, yp, zp, kind='linear'):
    interp = interpolate.interp2d(xp, np.log10(zp), np.log10(zp), kind=kind)
    def func(x, y): return np.power(10.0, interp(x, np.log10(y)))
    return func


def log_lin_log_interp_2d(xp, yp, zp, kind='linear'):
    interp = interpolate.interp2d(np.log10(xp), yp, np.log10(zp), kind=kind)
    def func(x, y): return np.power(10.0, interp(np.log10(x), y))
    return func


def log_log_lin_interp_2d(xp, yp, zp, kind='linear'):
    interp = interpolate.interp2d(np.log10(xp), np.log10(yp), zp, kind=kind)
    def func(x, y): return interp(np.log10(x), np.log10(y))
    return func


def log_log_log_interp_2d(xp, yp, zp, kind='linear'):
    interp = interpolate.interp2d(np.log10(xp), np.log10(yp), np.log10(zp), kind=kind)
    def func(x, y): return np.power(10.0, interp(np.log10(x), np.log10(y)))
    return func


def lin_lin_lin_interp_2d(xp, yp, zp, kind='linear'):
    interp = interpolate.interp2d(xp, yp, zp, kind=kind)
    def func(x, y): return interp(x, y)
    return func


# %% ELF utilities
def get_km_kp(E, hw):
    km = np.sqrt(2 * const.m / const.hbar ** 2) * (np.sqrt(E) - np.sqrt(E - hw))
    kp = np.sqrt(2 * const.m / const.hbar ** 2) * (np.sqrt(E) + np.sqrt(E - hw))
    return km, kp
