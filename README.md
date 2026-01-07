================================================================================
                           DEBRIS BROWSER FOR DEBIAN
================================================================================

Debris is a lightweight, GTK-based web browser designed for Debian-based Linux
distributions (Debian, Ubuntu, Mint, etc.). It uses WebKit2GTK for modern web
rendering and GStreamer for multimedia playback, providing a clean and minimal
browsing experience with low system resource usage.

--------------------------------------------------------------------------------
1. KEY FEATURES
--------------------------------------------------------------------------------

* MODERN ENGINE: Powered by WebKit2GTK for fast, standards-compliant browsing.
* GTK UI: Clean and simple interface built with GTK+ 3.
* MULTIMEDIA: GStreamer integration for audio and video playback.
* LIGHTWEIGHT: Minimal dependencies and fast startup.
* DEBIAN PACKAGE: Includes a .deb package for easy installation.

--------------------------------------------------------------------------------
2. PREREQUISITES (DEPENDENCIES)
--------------------------------------------------------------------------------

To compile Debris, install the required development libraries:

sudo apt update
sudo apt install g++ pkg-config libgtk-3-dev libwebkit2gtk-4.0-dev \
                 libgstreamer1.0-dev build-essential

--------------------------------------------------------------------------------
3. COMPILATION
--------------------------------------------------------------------------------

Use the following command to compile Debris from source:

g++ -o main main.cpp `pkg-config --cflags --libs gtk+-3.0 webkit2gtk-4.0 gstreamer-1.0`

This will generate the executable file "main".

--------------------------------------------------------------------------------
4. BUILDING THE DEBIAN PACKAGE
--------------------------------------------------------------------------------

To build the .deb package, run:

dpkg-deb --build Debris_0.5_deb

This will create the file:

Debris_0.5_deb.deb

--------------------------------------------------------------------------------
5. INSTALLATION
--------------------------------------------------------------------------------

Install the generated Debian package using:

sudo dpkg -i Debris_0.5_deb.deb

If any dependencies are missing, fix them with:

sudo apt --fix-broken install

--------------------------------------------------------------------------------
6. USAGE
--------------------------------------------------------------------------------

* Launch: Run "debris" from your terminal or application menu.
* Navigation: Standard browser controls (address bar, back/forward, reload).
* Multimedia: Supports HTML5 audio and video via GStreamer.

--------------------------------------------------------------------------------
7. TECHNICAL NOTES
--------------------------------------------------------------------------------

* Architecture: Built with GTK3, WebKit2GTK, and C++.
* OS Compatibility: Designed for Debian-based systems.
* Rendering: Uses WebKit2GTK for secure, multi-process browsing.
* Performance: Minimal UI overhead for fast startup and low memory usage.

--------------------------------------------------------------------------------
8. LICENSE
--------------------------------------------------------------------------------

Debris is free software: you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <https://www.gnu.org/licenses/>.

================================================================================

