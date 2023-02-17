#include <iostream>
#include <cmath>
#include <vector>
#include <random>
#include <ctime>
#include <SFML/Graphics.hpp>

// Constants
const int WINDOW_WIDTH = 800;
const int WINDOW_HEIGHT = 600;
const int PARTICLE_SIZE = 10;
const int NUM_PARTICLES = 50;
const float PARTICLE_SPEED = 100.0f;
const float PARTICLE_COLLISION_DAMPING = 0.8f;
const float WALL_COLLISION_DAMPING = 0.9f;
const int FPS = 60;

// Particle class
class Particle {
public:
    Particle(float x, float y, float size, float speed)
        : m_x(x), m_y(y), m_size(size), m_speed(speed), m_vx(0), m_vy(0)
    {
        // Set random velocity
        std::mt19937 rng(time(nullptr));
        std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
        m_vx = dist(rng) * m_speed;
        m_vy = dist(rng) * m_speed;
    }

    void update(float dt) {
        m_x += m_vx * dt;
        m_y += m_vy * dt;
    }

    void draw(sf::RenderWindow& window) {
        sf::CircleShape shape(m_size);
        shape.setFillColor(sf::Color::White);
        shape.setPosition(m_x - m_size, m_y - m_size);
        window.draw(shape);
    }

    void collideWithWall(int width, int height) {
        if (m_x - m_size < 0) {
            m_x = m_size;
            m_vx *= -WALL_COLLISION_DAMPING;
        } else if (m_x + m_size > width) {
            m_x = width - m_size;
            m_vx *= -WALL_COLLISION_DAMPING;
        }

        if (m_y - m_size < 0) {
            m_y = m_size;
            m_vy *= -WALL_COLLISION_DAMPING;
        } else if (m_y + m_size > height) {
            m_y = height - m_size;
            m_vy *= -WALL_COLLISION_DAMPING;
        }
    }

    void collideWithParticle(Particle& other) {
        float dx = other.m_x - m_x;
        float dy = other.m_y - m_y;
        float dist = std::sqrt(dx*dx + dy*dy);

        if (dist <= m_size + other.m_size) {
            // Calculate the angle of the collision
            float angle = std::atan2(dy, dx);

            // Rotate the velocities
            float self_vx_new = m_vx * std::cos(angle) + m_vy * std::sin(angle);
            float self_vy_new = m_vy * std::cos(angle) - m_vx * std::sin(angle);
            float other_vx_new = other.m_vx * std::cos(angle) + other.m_vy * std::sin(angle);
            float other_vy_new = other.m_vy * std::cos(angle) - other.m_vx * std::sin(angle);

            // Calculate the new velocities
            float self_vx_final = ((m_size - other.m_size) * self_vx_new + 2 * other.m_size * other_vx_new) / (m_size + other.m_size);
            float other_vx_final = ((other.m_size - m_size) * other_vx_new + 2 * m_size * self_vx_new) / (m_size + other
            float self_vy_final = self_vy_new;
            float other_vy_final = other_vy_new;

            // Rotate the velocities back
            m_vx = self_vx_final * std::cos(-angle) + self_vy_final * std::sin(-angle);
            m_vy = self_vy_final * std::cos(-angle) - self_vx_final * std::sin(-angle);
            other.m_vx = other_vx_final * std::cos(-angle) + other_vy_final * std::sin(-angle);
            other.m_vy = other_vy_final * std::cos(-angle) - other_vx_final * std::sin(-angle);

            // Move particles to avoid overlap
            float overlap = m_size + other.m_size - dist;
            m_x -= overlap/2 * std::cos(angle);
            m_y -= overlap/2 * std::sin(angle);
            other.m_x += overlap/2 * std::cos(angle);
            other.m_y += overlap/2 * std::sin(angle);

            // Dampen the velocities after collision
            m_vx *= PARTICLE_COLLISION_DAMPING;
            m_vy *= PARTICLE_COLLISION_DAMPING;
            other.m_vx *= PARTICLE_COLLISION_DAMPING;
            other.m_vy *= PARTICLE_COLLISION_DAMPING;
        }
    }

private:
    float m_x;
    float m_y;
    float m_size;
    float m_speed;
    float m_vx;
    float m_vy;
};

// Main function
int main() {
    // Create window
    sf::RenderWindow window(sf::VideoMode(WINDOW_WIDTH, WINDOW_HEIGHT), "Particle Simulation");
    window.setFramerateLimit(FPS);

    // Create particles
    std::vector<Particle> particles;
    for (int i = 0; i < NUM_PARTICLES; i++) {
        float x = std::rand() % (WINDOW_WIDTH - PARTICLE_SIZE*2) + PARTICLE_SIZE;
        float y = std::rand() % (WINDOW_HEIGHT - PARTICLE_SIZE*2) + PARTICLE_SIZE;
        particles.push_back(Particle(x, y, PARTICLE_SIZE, PARTICLE_SPEED));
    }

    // Start main loop
    sf::Clock clock;
    while (window.isOpen()) {
        // Handle events
        sf::Event event;
        while (window.pollEvent(event)) {
            if (event.type == sf::Event::Closed) {
                window.close();
            }
        }

        // Update particles
        float dt = clock.restart().asSeconds();
        for (int i = 0; i < NUM_PARTICLES; i++) {
            particles[i].update(dt);
            particles[i].collideWithWall(WINDOW_WIDTH, WINDOW_HEIGHT);

            for (int j = i+1; j < NUM_PARTICLES; j++) {
                particles[i].collideWithParticle(particles[j]);
            }
        }

        // Draw particles
        window.clear(sf::Color::Black);
        for (int i = 0; i < NUM_PARTICLES; i++) {
            particles[i].draw(window);
        }
        window.display();
    }

    return 0;
}
