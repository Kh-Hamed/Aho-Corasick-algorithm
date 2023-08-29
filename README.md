# Aho-Corasick-algorithm

Assume that we have a poem whose length is $100^{100}$. Try to find how beautiful is this poem. There are some good words(S<sub>i</sub>), with some weights assigned to each of them (C<sub>i</sub>). The formula for calculating the beauty of this poem string (T) is :

$$
\sum_{i=1}^{n} \frac{C_i \times \text{occurs}(T, S_i)}{|T|}
$$

Input: In the first line you have n, which shows the number of beautiful words. Then, in the next line, you'll have C<sub>i</sub> s. In the n following lines, you'll get n beautiful words.

<br>

|Sample Input | Output |
|:---------------------------------------------------|:-------:|
|1<br>10<br>aa                                       |10.000000|
|4<br>2<br>3<br>4<br>5<br>abb<br>bba<br>aab<br>baa   |3.500000|


**Note: Only around 16% of students could solve this problem 100% correctly for all of the test cases**

