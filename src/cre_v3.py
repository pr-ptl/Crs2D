#===============================================================================================================
#   2D COFFEE RING EFFECT
#===============================================================================================================

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from scipy.integrate import solve_ivp

#===============================================================================================================

class CRS:
    def __init__(self, R=1.0, h0=0.1, D=0.01, Pe=100, n_r=100, n_theta=100):
        """
        Initialize coffee ring simulator
        self : class pointer
        
        Parameters:
        - R: droplet radius
        - h0: initial droplet height at center
        - D: particle diffusion coefficient
        - Pe: Peclet number (advection/diffusion ratio)
        - n_r: number of radial grid points
        - n_theta: number of angular grid points
        """
        
        self.R = R
        self.h0 = h0
        self.D = D
        self.Pe = Pe
        self.n_r = n_r
        self.n_theta = n_theta
        
        # Spatial grid
        self.r = np.linspace(0.01, R, n_r)  # Exempting singularity at r=0
        self.theta = np.linspace(0, 2*np.pi, n_theta)
        self.dr = self.r[1] - self.r[0]
        
        # Droplet profile (spherical cap)
        self.h = self.cap_prf(self.r)
        
        # Particle concentration (uniform)
        self.c = np.ones_like(self.r)
        
        # Initialize deposited particles
        self.deposited = np.zeros_like(self.r)
        
        # Time parameters
        self.t = 0
        self.dt = 0.001
        
    def cap_prf(self, r):
        
        """Calculate spherical cap height profile"""
        
        # Contact angle
        ca = np.pi/3
        
        # Profile height
        Rc = self.R / np.sin(ca)
        h_edge = self.R * np.tan(np.pi/2 - ca)
        
        disc = Rc**2 - r**2
        #disc = np.maximum(disc, 0.01)
        
        h = np.sqrt(disc) - (Rc - h_edge)
        #h = np.maximum(h, 0.01)
        
        return h
    
    def evap_ratef(self, r):
        
        """Evaporation rate as function of radius"""
        
        ee = 1 + 10 * np.exp(-10 * (self.R - r))
        base_rate = 0.05
        
        return base_rate * ee
    
    def vr(self, r):
        
        """Calculate radial velocity due to evaporation-driven flow"""
        
        if r <= 0.01:
            return 0.0
            
        r_int = np.linspace(0.01, r, 100)  
        evap_vals = [self.evap_ratef(ri) for ri in r_int]
        Jint = np.trapz(evap_vals, r_int)
        
        # Velocity from continuity: v_r = J_integral / h
        
        idx = np.argmin(np.abs(self.r - r))
        v_r = Jint / self.h[idx]
        
        return v_r
    
    def update_height(self):
        
        """Update droplet height due to evaporation"""
        
        evap_rate = np.array([self.evap_ratef(ri) for ri in self.r])
        self.h -= evap_rate * self.dt
        self.h = np.maximum(self.h, 0.01)
    
    def update_conc(self):
        
        """Update particle concentration using advection-diffusion equation"""
        
        c_new = self.c.copy()
        
        for i in range(1, len(self.r)-1):
            r_i = self.r[i]
            h_i = self.h[i]
            
            # Calculate radial velocity
            v_r = self.vr(r_i)
            
            # Advection term: -v_r * dc/dr
            dc_dr = (self.c[i+1] - self.c[i-1]) / (2 * self.dr)
            adv = -v_r * dc_dr
            
            # Diffusion term: D * (1/r) * d/dr(r * dc/dr)
            dc_dr_f = (self.c[i+1] - self.c[i]) / self.dr
            dc_dr_b = (self.c[i] - self.c[i-1]) / self.dr
            
            r_f = self.r[i] + self.dr/2
            r_b = self.r[i] - self.dr/2
            
            diff = self.D * (r_f * dc_dr_f - r_b * dc_dr_b) / (r_i * self.dr)
            
            # Concentration due to evaporation (conservation)
            evap_conc = self.c[i] * self.evap_ratef(r_i) / h_i
            
            # Update concentration
            c_new[i] += self.dt * (diff + adv + evap_conc)
        
        # Boundary conditions
        c_new[0] = c_new[1]  
        c_new[-1] = max(0, c_new[-1])
        
        self.c = np.maximum(c_new, 0)
    
    def update_depo(self):
        
        """Update deposited particles"""
        
        depo_rate = 0.05 * self.c * np.exp(-self.h)
        
        edge_factor = 1 + 5 * np.exp(-5 * (self.R - self.r))
        depo_rate *= edge_factor
        
        depo_amount = depo_rate * self.dt
        self.deposited += depo_amount
        self.c -= depo_amount / np.maximum(self.h, 0.01)
        self.c = np.maximum(self.c, 0)
    
    def step(self):
        
        """Perform one time step of the simulation"""
        
        self.update_height()
        self.update_conc()
        self.update_depo()
        self.t += self.dt
    
    def simulate(self, total_time=10.0, save_interval=50):
        
        times = []
        heights = []
        conc = []
        depo = []
        
        n_steps = int(total_time / self.dt)
        
        for i in range(n_steps):
            if i % save_interval == 0:  # Store data
                times.append(self.t)
                heights.append(self.h.copy())
                conc.append(self.c.copy())
                depo.append(self.deposited.copy())
                
                if i % (save_interval * 10) == 0:
                    print(f"Progress: {i/n_steps:.1f}%, t={self.t:.3f}")
            
            self.step()
            
            if np.max(self.h) < 0.02:
                print(f"Droplet evaporated at t={self.t:.3f}")
                break
        
        return np.array(times), np.array(heights), np.array(conc), np.array(depo)

def plot_fdepo(sim, depo):

    fig, ax = plt.subplots(figsize=(10, 6))
    ax.plot(sim.r/sim.R, depo[-1], 'b-', linewidth=3, label='Final Deposition')
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'$\sigma$', fontsize=14)
    ax.set_title('Coffee Ring Pattern', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=12)
    plt.tight_layout()
    
    return fig

def plot_hevol(sim, times, heights):
    
    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(0, len(times), max(1, len(times)//6)):
        ax.plot(sim.r/sim.R, heights[i], alpha=0.8, linewidth=2,
                label=f't={times[i]:.2f}')
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'Height $h$', fontsize=14)
    ax.set_title('Evolution of Droplet Height', fontsize=16, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig

def plot_concevol(sim, times, conc):
    
    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(0, len(times), max(1, len(times)//6)):
        ax.plot(sim.r/sim.R, conc[i], alpha=0.8, linewidth=2,
                label=f't={times[i]:.2f}')
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'Concentration $C$', fontsize=14)
    ax.set_title('Evolution of Particle Concentration', fontsize=16, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    
    return fig

def plot_depoevol(sim, times, depo):
    
    fig, ax = plt.subplots(figsize=(10, 6))
    for i in range(0, len(times), max(1, len(times)//6)):
        ax.plot(sim.r/sim.R, depo[i], alpha=0.8, linewidth=2,
                label=f't={times[i]:.2f}')
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'$\sigma$', fontsize=14)
    ax.set_title('Evolution of Particle Deposition', fontsize=16, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    return fig

def plot_2df(sim, deposited):
    
    x = np.linspace(-sim.R, sim.R, 200)
    y = np.linspace(-sim.R, sim.R, 200)
    X, Y = np.meshgrid(x, y)
    R_grid = np.sqrt(X**2 + Y**2)
    
    Z = np.zeros_like(R_grid)
    for i, r_val in enumerate(sim.r):
        mask = (R_grid >= r_val - sim.dr/2) & (R_grid < r_val + sim.dr/2)
        Z[mask] = deposited[i]
    
    Z[R_grid > sim.R] = 0
    
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(Z, extent=[-sim.R, sim.R, -sim.R, sim.R], 
                   origin='lower', cmap='YlOrBr')
    ax.set_title('Coffee Ring', fontsize=16, fontweight='bold')
    ax.set_xlabel(r'$x/R$', fontsize=14)
    ax.set_ylabel(r'$y/R$', fontsize=14)
    cbar = plt.colorbar(im, label=r'$\sigma$')
    cbar.ax.tick_params(labelsize=12)
    
    circle = plt.Circle((0, 0), sim.R, fill=False, color='black', linewidth=3)
    ax.add_patch(circle)
    
    plt.tight_layout()
    
    return fig

def animate_hevol(sim, times, heights):
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    line, = ax.plot([], [], 'b-', linewidth=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, np.max(heights) * 1.1)
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'Height $h$', fontsize=14)
    ax.set_title('Evolution of Droplet Height', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=14,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    def animate(frame):
        if frame >= len(times):
            frame = len(times) - 1
        line.set_data(sim.r/sim.R, heights[frame])
        time_text.set_text(f'Time: {times[frame]:.3f}')
        return line, time_text
    
    anim = FuncAnimation(fig, animate, frames=len(times), interval=100, 
                        blit=True, repeat=True)
    plt.tight_layout()
    
    return fig, anim

def animate_concevol(sim, times, conc):
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    line, = ax.plot([], [], 'r-', linewidth=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, np.max(conc) * 1.1)
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'Concentration $C$', fontsize=14)
    ax.set_title('Evolution of Particle Concentration', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=14,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    def animate(frame):
        if frame >= len(times):
            frame = len(times) - 1
        line.set_data(sim.r/sim.R, conc[frame])
        time_text.set_text(f'Time: {times[frame]:.3f}')
        return line, time_text
    
    anim = FuncAnimation(fig, animate, frames=len(times), interval=100, 
                        blit=True, repeat=True)
    plt.tight_layout()
    
    return fig, anim

def animate_depoevol(sim, times, depo):
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    line, = ax.plot([], [], 'g-', linewidth=3)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, np.max(depo[-1]) * 1.1)
    ax.set_xlabel(r'$r/R$', fontsize=14)
    ax.set_ylabel(r'$\sigma$', fontsize=14)
    ax.set_title('Evolution of Particle Deposition', fontsize=16, fontweight='bold')
    ax.grid(True, alpha=0.3)
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=14,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    def animate(frame):
        if frame >= len(times):
            frame = len(times) - 1
        line.set_data(sim.r/sim.R, depo[frame])
        time_text.set_text(f'Time: {times[frame]:.3f}')
        return line, time_text
    
    anim = FuncAnimation(fig, animate, frames=len(times), interval=100, 
                        blit=True, repeat=True)
    plt.tight_layout()
    
    return fig, anim

def animate_2devol(sim, times, depo):
    
    fig, ax = plt.subplots(figsize=(10, 6))

    x = np.linspace(-sim.R, sim.R, 200)
    y = np.linspace(-sim.R, sim.R, 200)
    X, Y = np.meshgrid(x, y)
    R_grid = np.sqrt(X**2 + Y**2)

    Z_init = np.zeros_like(R_grid)
    im = ax.imshow(Z_init, extent=[-sim.R, sim.R, -sim.R, sim.R], 
                   origin='lower', cmap='YlOrBr', vmin=0, vmax=np.max(depo[-1]))
    
    ax.set_title('Coffee Ring', fontsize=16, fontweight='bold')
    ax.set_xlabel(r'$x/R$', fontsize=14)
    ax.set_ylabel(r'$y/R$', fontsize=14)
    
    cbar = plt.colorbar(im, label=r'$\sigma$')
    cbar.ax.tick_params(labelsize=12)
    
    circle = plt.Circle((0, 0), sim.R, fill=False, color='white', linewidth=3)
    ax.add_patch(circle)
    
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes, fontsize=14,
                       bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    
    def animate(frame):
        if frame >= len(times):
            frame = len(times) - 1
            
        Z = np.zeros_like(R_grid)
        for i, r_val in enumerate(sim.r):
            mask = (R_grid >= r_val - sim.dr/2) & (R_grid < r_val + sim.dr/2)
            Z[mask] = depo[frame][i]
        Z[R_grid > sim.R] = 0
        
        im.set_data(Z)
        time_text.set_text(f'Time: {times[frame]:.3f}')
        
        return [im, time_text]
    
    anim = FuncAnimation(fig, animate, frames=len(times), interval=150, 
                        blit=False, repeat=True)
    plt.tight_layout()
    return fig, anim

if __name__ == "__main__":
    
    sim = CRS(R=1.0, h0=0.1, Pe=100, n_r=100)
    
    print("Running coffee ring simulation...")
    print("This models evaporation-driven flow and particle deposition")
    
    times, heights, conc, depo = sim.simulate(total_time=10.0, save_interval=20)
    
    print(f"Simulation completed after {len(times)} time steps")
    print(f"Final time: {times[-1]:.3f}")
    
    print("Creating plots...")
    fig1 = plot_fdepo(sim, depo)
    fig2 = plot_hevol(sim, times, heights)
    fig3 = plot_concevol(sim, times, conc)
    fig4 = plot_depoevol(sim, times, depo)
    fig5 = plot_2df(sim, depo[-1])
    
    print("Creating animations...")
    fig6, anim1 = animate_hevol(sim, times, heights)
    fig7, anim2 = animate_concevol(sim, times, conc)
    fig8, anim3 = animate_depoevol(sim, times, depo)
    fig9, anim4 = animate_2devol(sim, times, depo)
    
    plt.show()
    
    save_animations = input("Save animations as GIF files? (y/n): ").lower() == 'y'
    if save_animations:
        print("Saving animations...")
        anim1.save('height_evol.gif', writer='pillow', fps=8)
        anim2.save('conc_evol.gif', writer='pillow', fps=8)
        anim3.save('depo_evol.gif', writer='pillow', fps=8)
        anim4.save('2d_ring.gif', writer='pillow', fps=6)
        print("Animations saved!")
    
    # Statistics
    fdepo = depo[-1]
    cdepo = np.mean(fdepo[:len(fdepo)//4])
    edge_depo = np.mean(fdepo[-len(fdepo)//4:])
    
    print(f"\nCoffee Ring Statistics:")
    print(f"Center deposition: {cdepo:.3f}")
    print(f"Edge deposition: {edge_depo:.3f}")
    print(f"Ring enhancement factor: {edge_depo/cdepo:.2f}")
