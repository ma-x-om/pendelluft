## The Problem
You might be wondering "why so many lists in Code_3.py and what's the positioning logic behind them?". So let me explain myself:

In Victor's original MATLAB codes, every single parameter was called out in a separated variable such as LEmin, REmax, L1min, L2max and so goes...
It was an eye sore. It needed to be changed. We had a perfect set of images and potentials and everything else and we could well nest those parameters in an array or a list.

<i>In the original code:</i><br>
![image](https://user-images.githubusercontent.com/126175949/223301285-e54f5dc7-899c-4f0c-8414-a8c2233b76f8.png)

## Disambiguation
Since lists can be called out as empty lists and we could easily append items to them we chose to use them over NumPy arrays. Our next mission was to understand how to nest together those parameters. We needed to divide left from right and separate potE (Divergence-free) and potW (Curl-free) potentials.

<i> Used settings to distingish between Left, Right, E and W</i><br>
![image](https://user-images.githubusercontent.com/126175949/223301936-02e11647-65ce-4df5-9cb4-b7b974a30812.png)
![image](https://user-images.githubusercontent.com/126175949/223302005-95682b08-831f-4be6-8950-81577d8d91ef.png)
![image](https://user-images.githubusercontent.com/126175949/223302110-6d2ae112-3977-4c82-904e-c6d4c1adc89d.png)

## The Solution
The chosen method was to nest parameters by TYPE and organize by QUADRANT. As shown in the previous pictures, the QUADRANTS organization follows as Q1 means LEFT-MAX, Q2 RIGHT-MAX, Q3 LEFT-MIN and Q4 RIGHT-MIN:

| QUADRANT | LEFT | RIGHT |
| :------: | :--: | :---: |
| MAX | 1 | 2 |
| MIN | 3 | 4 |

And the parameters are stored in 6 differents lists/arrays by type:
- arr_LE = [LE1, LE2, LE3, LE4]
- arr_LW = [LW1, LW2, LW3, LW4]
- arr_E = [E1, E2, E3, E4]
- arr_W = [W1, W2, W3, W4]
- arr_LEm = [LEm1, LEm2, LEm3, LEm4]
- arr_LWm = [LWm1, LWm2, LWm3, LWm4]

That way we can aliviate problemas related to array dimensions and iterate over the lists and arrays with a relative ease compared to before.

### Other possible solution
![image](https://user-images.githubusercontent.com/126175949/223304578-6d2338c2-033f-463a-9810-725eca998f46.png)
