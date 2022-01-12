# LangrangianDoublePendulum

**[CLICK HERE TO SEE IT IN ACTION](https://tailcalled.github.io/LangrangianDoublePendulum/double_pendulum.html)**

A double pendulum simulation I made using Lagrangian mechanics in Javascript. It's kinda broken; probably partly this is due to numerical precision errors in the differentiation, equation solving and integration, but it might also just be implemented wrong to a degree.

Basic principles:

 * A Langrangian for objects with mass in a gravitational field is computed.
 * This is transformed into a Langrangian for a double pendulum, by using a function that maps angles to x/y positions, and transforming the velocities using numerical differentation and the chain rule.
 * The Euler-Lagrange equation is set up, again using numerical differentiation.
 * Using Newton's method, numerical differentation and linear algebra, we solve the Euler-Lagrange equation.
 * This gets integrated and rendered, with some extra debug info.

The code is completely standalone and uses no libraries. Warning, it is kind of ugly because it's just something I quickly threw together.
