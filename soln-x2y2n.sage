# The below function returns the minimum of 't' unrelated solutions (i.e., x^2+y^2=n and 0<=x<=y) or all the unrelated solutions, by doing a brute force search on x. 

def nontrivsolnx2y2n(n,t):
	lst = []
	count = 0
	x = 0
	sqrtnb2 = floor(sqrt(n/2))
	while count <= t and x <= sqrtnb2:
		tmp1 = sqrt((n-x^2))
		tmp2 = Integer(floor(tmp1))
		if tmp1 == tmp2:
			lst.append((x,tmp2))
			count += 1
		x += 1
	return lst


# The below function implements the reduction of factoring of n (n >= 3) to the problem enumerating at most 3 unrelated solutions to x^2+y^2=n. 

def factorenum(n):
	if n%2 == 0:
		print "Non-trivial factor of n is ", 2
		return
	if Integer(n).is_prime(proof=True):
		print "n is a prime."
		return
	tmp = Integer(n).perfect_power()
	if tmp[0] < n:
		print "n is a perfect power. The base and exponent are ", tmp
		return
	lst = nontrivsolnx2y2n(n,2)	# Proof reduction needs t=3, but t=2 works in practice.
	if len(lst) == 0:
		print "There is no solution for x^2+y^2=n. Hence cannot factor n."
		return
	for i in range(len(lst)):
		tmpx = gcd(lst[i][0],n)
		tmpy = gcd(lst[i][1],n)		
		if tmpx > 1:
			print "Non-trivial factor of n via gcd with x is ", tmpx
			return
		elif tmpy > 1:
			print "Non-trivial factor of n via gcd with y is ", tmpy
			return
	zn = IntegerModRing(n)
	tmpr1 = Integer(zn(lst[0][0])/zn(lst[0][1]))
	for i in range(1,len(lst)):
		tmpr2 = Integer(zn(lst[i][0])/zn(lst[i][1]))
		if tmpr2 != tmpr1 and tmpr2 != n-tmpr1:
			tmp = gcd(tmpr1-tmpr2,n)
			if tmp > 1:
				print "Non-trivial factor of n via differnce of roots is ", tmp
				return
			else:
				print "Error: gcd=1"	
	print "Error: all cases covered"	
	return		


##### Auxilliary routines for investigating various properties of solutions to x^2+y^2=n. ######

def genpmod4_1(n):
	tlst = prime_range(3,n)
	lst = []
	for i in range(len(tlst)):
		if tlst[i]%4 == 1:
			lst.append(tlst[i])
	return lst

def comproot(a,e,n):
	lst = []
	for i in range(n):
		if (i^e)%n == a%n:
			lst.append(i)
	return lst


def compsolnroot_1(n):
	zn = IntegerModRing(n)
	tmplst1 = comproot(-1,2,n)
	print "Sq roots of -1: ", tmplst1
	tmplst2 = nontrivsolnx2y2n(n,n)
	print "Non trivial solns: ", tmplst2
	lst1 = []
	lst2 = []
	for i in tmplst2:
		if gcd(i[0],n)>1 or gcd(i[1],n)>1:
			lst1.append(0)		# Non-trivial factor of n has been/can be found
			lst2.append(n)
		else:
			a = Integer(zn(i[0])/zn(i[1]))
			lst1.append(a)
			lst2.append(n-a)
	print "Sq. roots of -1 (mod n) and their complement via xy^-1"
	print lst1
	print lst2
	 	

def enumcompsolnroot_1(n):
	plst = genpmod4_1(n)
	for i in range(len(plst)-1):
		for j in range(i+1,len(plst)):
			p = plst[i]*plst[j]
			compsolnroot_1(p)
		

