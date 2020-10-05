There are two files to validate our work:
1) attack.py
2) soln-x2y2n.sage


1) The file attack.py consists of our rider's attack proposed in Section 2.2 and Algorithm 1.

The code (written in Python) uses Google Maps API and the UTM library for python.
Due to reasons of confidentiality, we are unable to provide our API key along with the code.

The API key is necessary for using Google Map API's services. The API key is linked with a billing account,
and can be set up by following the instructions here:
https://developers.google.com/maps/documentation/javascript/get-api-key

For our code, we use the Road API and Distance-Matrix API, and these need to be activated from the 
Google Cloud Platform Console. 

2) The file soln-x2y2n.sage computationally verifies the proof of Theorem 4.2. 
More specifically, it enumerates two non-trivial solutions (x1, y1), (x2, y2) to x^2+y^2=N and verifies
that x1y2 + x2y1 < N
