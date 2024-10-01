#include <iostream>
#include <stdexcept>
#include <limits.h>
#include <cmath>

const int SIGN_BIT_POSITION = 31;
const int B_POSITION = 26;
const int MAX_BITS = 26;
const int MAX_B_VALUE = 31;
const int ARRAY_SIZE = 1000000;     
const double TOLERANCE = 1e-5;

int gcd(int a, int b) {
    if (a < b) {
        return gcd(b, a);
    }
    else if (b != 0) {
        return gcd(b, a % b);
    }
    else {
        return a;
    }
}

int numBits(int n) {
    int res = 0;
    while (n != 0) {
        n >>= 1; ++res;
    }
    return res;
}

void printBits(unsigned num) {

    bool negative = num & (1 << SIGN_BIT_POSITION);
    int B = (num >> B_POSITION) % (MAX_B_VALUE + 1);
    int numerator = (num >> B) % (1 << (MAX_BITS - B));
    int denominator = num % (1 << B);

    std::cout << negative << " ";

    for (int i = 4; i >= 0; --i) {
        std::cout << ((B >> i) & 1);
    }
    std::cout << " ";

    for (int i = 25; i >= B; --i) {
        std::cout << ((numerator >> (i - B)) & 1);
    }
    std::cout << " ";

    for (int i = B - 1; i >= 0; --i) {
        std::cout << ((denominator >> i) & 1);
    }
    std::cout << std::endl;
}




unsigned makeRat(int n, int d) {
    if (d == 0) {
        throw std::invalid_argument("Denominator cannot be zero.");
    }

    bool negative = (n < 0) ^ (d < 0);
    n = abs(n);
    d = abs(d);

    int gcdND = gcd(n, d);
    if (gcdND != 1)
    {
        n /= gcdND;
        d /= gcdND;
    }
    
    int B = numBits(d);

    if (B == MAX_B_VALUE || (n >> (MAX_BITS - B)) > 0) {
        throw std::overflow_error("Numerator or denominator too large to fit in 26 bits.");
    }

    unsigned result = 0;
    if (negative) {
        result |= (1 << SIGN_BIT_POSITION);
    }
    result |= (B << B_POSITION);     
    result |= d % (1 << B);       
    result |= (n << B);      

    return result;
}

int getDenom(unsigned rat) {

    int B = (rat >> B_POSITION) % (MAX_B_VALUE + 1); 

    int denominator = rat % (1 << B);

    return denominator;
}

//BONUS
/*
unsigned makeRat(int n, int d) {
    if (d == 0) {
        throw std::invalid_argument("Denominator cannot be zero.");
    }

    bool negative = (n < 0) ^ (d < 0);
    n = abs(n);
    d = abs(d);

    int gcdND = gcd(n, d);
    if (gcdND != 1) {
        n /= gcdND;
        d /= gcdND;
    }

    int B = numBits(d) - 1;

    if (B == MAX_B_VALUE || (n >> (MAX_BITS - B)) > 0) {
        throw std::overflow_error("Numerator or denominator too large to fit in 26 bits.");
    }

    unsigned result = 0;
    if (negative) {
        result |= (1 << SIGN_BIT_POSITION);
    }
    result |= (B << B_POSITION);
    result |= d % (1 << B);
    result |= (n << B);

    return result;
}

int getDenom(unsigned rat) {
    int B = (rat >> B_POSITION) % (MAX_B_VALUE + 1);

    int denominator = rat % (1 << B);

    denominator |= (1 << B);

    return denominator;
}
*/
int getNumer(unsigned rat) {
    bool negative = rat & (1 << SIGN_BIT_POSITION);

    int B = (rat >> B_POSITION) % (MAX_B_VALUE + 1);

    int numerator = (rat >> B) % (1 << (MAX_BITS - B));

    if (negative) {
        numerator = -numerator;
    }

    return numerator;
}

unsigned ratAdd(unsigned r1, unsigned r2) {
    int n1 = getNumer(r1);
    int d1 = getDenom(r1);
    int n2 = getNumer(r2);
    int d2 = getDenom(r2);

    int numerator = n1 * d2 + n2 * d1;
    int denominator = d1 * d2;

    return makeRat(numerator, denominator);
}

unsigned ratNeg(unsigned rat) {
    return rat ^ (1 << SIGN_BIT_POSITION);
}

unsigned ratAbs(unsigned rat) {
    return rat & ~(1 << SIGN_BIT_POSITION);
}

unsigned ratSub(unsigned r1, unsigned r2) {
    return ratAdd(r1, ratNeg(r2));
}

unsigned ratMul(unsigned r1, unsigned r2) {
    int n1 = getNumer(r1);
    int d1 = getDenom(r1);
    int n2 = getNumer(r2);
    int d2 = getDenom(r2);

    int numerator = n1 * n2;
    int denominator = d1 * d2;

    return makeRat(numerator, denominator);
}

unsigned ratRec(unsigned rat) {
    int numerator = getNumer(rat);
    if (numerator == 0) {
        throw std::exception("Cannot reciprocate a zero numerator.");
    }

    bool negative = rat & (1 << SIGN_BIT_POSITION);

    unsigned denominator = getDenom(rat);

    unsigned newN = denominator;
    unsigned newD = numerator;

    return makeRat(newN, newD);
}


unsigned ratDiv(unsigned r1, unsigned r2) {
    return ratMul(r1, ratRec(r2));
}

float rat2float(unsigned rat) {
    int numerator = getNumer(rat);
    int denominator = getDenom(rat);

    float result = static_cast<float>(numerator) / static_cast<float>(denominator);

    return result;
}

void compareSums(unsigned rats[ARRAY_SIZE]) {
    double floatingBarSum = 0.0;
    float fp32Sum = 0.0f; 

    for (int i = 0; i < 4; i++) {
        int numerator = getNumer(rats[i]);
        int denominator = getDenom(rats[i]);

        floatingBarSum += static_cast<double>(numerator) / static_cast<double>(denominator);
        fp32Sum += rat2float(rats[i]);
    }

    std::cout << "floating bar sum: " << floatingBarSum << std::endl;
    std::cout << "fp32 sum: " << fp32Sum << std::endl;

    double difference = std::abs(floatingBarSum - fp32Sum);
    std::cout << "Difference: " << difference << std::endl;

    std::cout << (difference < TOLERANCE ? "The sums are very close to each other."
                                         : "There is a significant difference between the sums.")
        << std::endl;
}

int main() {
    try {
        unsigned r = makeRat(1, 3);
        unsigned r1 = makeRat(9, 8);
        unsigned r2 = makeRat(-3, 2);
        unsigned r3 = makeRat(0, 1);        
                
        std::cout << getNumer(r) << "/" << getDenom(r) << " = " << rat2float(r) << " : "; printBits(r);
        std::cout << getNumer(r1) << "/" << getDenom(r1) << " = " << rat2float(r1) << " : "; printBits(r1);
        std::cout << getNumer(r2) << "/" << getDenom(r2) << " = " << rat2float(r2) << " : "; printBits(r2);
        std::cout << getNumer(r3) << "/" << getDenom(r3) << " = " << rat2float(r3) << " : "; printBits(r3);
        std::cout << "Reciprocal numbers:" << std::endl;
        std::cout << getNumer(ratRec(r)) << "/" << getDenom(ratRec(r)) << std::endl;
        std::cout << getNumer(ratRec(r1)) << "/" << getDenom(ratRec(r1)) << std::endl;
        std::cout << getNumer(ratRec(r2)) << "/" << getDenom(ratRec(r2)) << std::endl;
        //std::cout << getNumer(ratRec(r3)) << "/" << getDenom(ratRec(r3)) << std::endl; // throw exeption
        std::cout << "Aritmetical operations:" << std::endl;
        std::cout << getNumer(r) << "/" << getDenom(r) << " + " << getNumer(r1) << "/" << getDenom(r1) << " = " << getNumer(ratAdd(r, r1)) << "/" << getDenom(ratAdd(r, r1)) << std::endl;
        std::cout << getNumer(r) << "/" << getDenom(r) << " - " << getNumer(r1) << "/" << getDenom(r1) << " = " << getNumer(ratSub(r, r1)) << "/" << getDenom(ratSub(r, r1)) << std::endl;
        std::cout << getNumer(r) << "/" << getDenom(r) << " * " << getNumer(r1) << "/" << getDenom(r1) << " = " << getNumer(ratMul(r, r1)) << "/" << getDenom(ratMul(r, r1)) << std::endl;
        std::cout << getNumer(r) << "/" << getDenom(r) << " / " << getNumer(r1) << "/" << getDenom(r1) << " = " << getNumer(ratDiv(r, r1)) << "/" << getDenom(ratDiv(r, r1)) << std::endl;
        std::cout << "Abs of " << getNumer(r2) << "/" << getDenom(r2) << " = " << getNumer(ratAbs(r2)) << "/" <<getDenom(ratAbs(r2)) << std::endl;
        std::cout << "Neg of " << getNumer(r) << "/" << getDenom(r) << " = " << getNumer(ratNeg(r)) << "/" <<getDenom(ratNeg(r)) << std::endl;

    }
    catch (const std::overflow_error& e) {
        std::cerr << "Overflow error: " << e.what() << std::endl;
    }
    catch (const std::invalid_argument& e) {
        std::cerr << "Invalid argument: " << e.what() << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "An unexpected error occurred: " << e.what() << std::endl;
    }

    return 0;
}




//Напишете в кратък коментар в кода по едно предимство / недостатък на това представяне.
// Предимство:
// -С floating bar можем да представим точно рационални числа, 
//  без закръгляне или загуба на точност, което е често срещано при представянето с float или double.
// Недостатък:
// -Ограничение в размерите на числителя и знаменателя: 
//  методът има ограничение в общия брой на битовете за числителя и знаменателя (26 бита),
//  което може да бъде недостатъчно за представяне на дроби с големи числители или знаменатели.
// Пример: 1/3 с float се представя като 0.33333333…, което е приближено и не е точно.
// Пример: 1200000000 няма да се побере в 26 бита
//Напишете по едно предимство/недостатък на това усложнение.
// Предимство:
// -По-ефективно използване на битовете: 
//  Тъй като знаменателят винаги има водещ бит единица, това позволява оптимизация в броя на битовете, 
//  необходими за представянето на знаменателя. Така се освобождават повече битове за числителя,
//  което позволява по-големи дроби да бъдат представени с дадения брой битове.
// Недостиг:
// -усложнява процеса на декодиране на числото.
//  сега е необходимо да се добави допълнителна логика за възстановяване на пълния знаменател, 
//  което увеличава сложността на функциите за работа с такива числа и може да доведе до по-високи изчислителни разходи.
