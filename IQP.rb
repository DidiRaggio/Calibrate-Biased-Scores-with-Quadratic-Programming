require 'simplex'
require 'matrix' #requiere Matrix

#generate radom dumber within range
def range (min, max)
    rand * (max-min) + min
end

#Generate array with with specified elements and fill out length with zeros
def fixed_array(size, other)  
   Array.new(size) do |i| 
   	if other[i] != nil
   		other[i]
   	else
   		0
   	end
   end
end

#generate zeros array of specified length
def empty_array(size)  
   Array.new(size) { |i| 0 }
end

#extend matrix class
class Matrix
  def push_to_matrix(i, j, x)
    @rows[i-1][j-1] = x
  end
end

#checks that all values are equal to ture
class Array
  def same_true_values?
    self.uniq.length == 1 && self[0] == true
    
  end
end

#modify "a" just ot include the columns specified in S for current iteration
def current_a(s, a0)
	a = Matrix[]
	s.each do |n|
		a =  Matrix.rows(a.to_a << a0.row(n-1))
	end
	a
end

#modify "b" just ot include the columns specified in S for current iteration
def current_b(s, b0)
	b = Vector[]
	s.each do |n|
		b =  Vector.elements(b.to_a << b0.component(n-1))
	end
	b
end

#generate KKT matrix
def generate_KKT_matrix(p, a)
	if a.column_count != 0
		k1 = Matrix[]
		p.row(0).size.times do |n|
			k1 =  Matrix.rows(k1.to_a << (p.row(n).to_a.concat a.transpose.row(n).to_a))
		end
		k2 = Matrix[]
		a.column(0).size.times do |n|
			k2 =  Matrix.rows(k2.to_a << (a.row(n).to_a.concat empty_array(a.column(0).size)))
		end
		k = k1.vstack(k2)
	else
		k = p
	end
end


def check_feasability(a0, xs, b0)
	feasability_check = a0*xs
	feasable = (feasability_check.to_a.zip(b0.to_a).map { |x, y| x <= y }).same_true_values?
	feasable
end

#check if x star is within the feasable space for the constraints
# def check_feasability(a0, xs, b0)
# 	feasability_check = a0*xs
# 	feasable = nil
# 	b0.size.times do |n|
# 		if feasability_check[n] <= b0[n] 
# 			feasable = true
# 		else
# 			feasable = false
# 		end
# 	end
# 	feasable
# end

#check if x star is an optimizer
def check_optimizer(nus)
	
	optimal = (nus.to_a.map { |x| x >= 0 }).same_true_values?
	optimal
end

#initialize matrixes and vectors
p = Matrix[[2, 0], [0, 2]]
q = Vector[-8, -6]
a = Matrix[[-1, 0], [0, -1], [1, 1]]
b = Vector[0, 0, 5]


# p = Matrix[[1, 0], [0, 3/2.0]]
# q = Vector[-1, -1]
# a = Matrix[[-1, -1], [-1, 0], [0, -1]]
# b = Vector[-6, 0, 0]


iteration_number = 0


#set initial conditions
a0 = a.clone
b0 = b.clone
s = Vector.elements([*1..(a0.row_count-1)])
s0 = s.clone
nus = Vector[]
x0 = current_a(s, a0).lup.solve(current_b(s, b0))#.elements(empty_array(a0.row_count-1))
prev_x = x0.clone
xs = x0.clone




while !check_optimizer(nus)

	while check_feasability(a0, xs, b0) && !check_optimizer(nus)
		#Iteration according to Karush–Kuhn–Tucker conditions 

		#update data for iteration
		iteration_number += 1
		prev_x = xs.clone

		#set current iteration matrixes
		a = current_a(s, a0)
		b = current_b(s, b0)
		k = generate_KKT_matrix(p, a)

		# concat qb
		qb = (q*-1).to_a.concat b.to_a
			# l, u, p = k.lup
			# puts l.lower_triangular? # => true
			# puts u.upper_triangular? # => true
			# puts p.permutation?      # => true
			# puts l * u == p * k      # => true

		#solve matrix system
		result = k.lup.solve(qb)

		#x star EQP
		xs = Vector.elements(result[0..(x0.size-1)])
		#nu star EQP
		nus = Vector.elements(result[(x0.size)..result.size])

		# set new constratints set for next iteration
		if s.size != 0 
			s_next = []
			s_next += s.to_a
			search_index = nus.to_a.bsearch_index{|x| x < 0}
			if search_index != nil
				s_next.delete_at(search_index)
				s =  Vector.elements(s_next)
			else
			s =  Vector[]
			end
		end

		puts "ITERATION #{iteration_number}"
		puts "The value of s: #{s}"
		puts "The value of a: #{a}"
		puts "The value of k: #{k}"
		puts "The value of q: #{q}"
		puts "The value of b: #{b}"
		puts "The value of resutl: #{result}"
		puts "The value of xs: #{xs}"
		puts "The value of nus: #{nus}"
		puts "The value of prev_x: #{prev_x}"
		puts "END OF ITERATION"

		puts "Current status: is feasable: #{check_feasability(a0, xs, b0)}, is optimal: #{check_optimizer(nus)}"


	end

	if !check_optimizer(nus)
		# Maximize the t function and modify constraints
		a0_prime = Matrix[]
		a0_beta = (a0 * ( xs - prev_x )).map(&:to_i)
		a0_beta.size.times do |n|
			a0_prime =  Matrix.rows(a0_prime.to_a << [a0_beta.component(n)])
		end
		b0_prime = (b0 - ( a0 * prev_x )).map(&:to_i)

		#solves simplex equation and obtains the max value of "t"
		simplex = Simplex.new(
			[1],
			a0_prime.to_a,
			b0_prime.to_a,
		)

		s = Vector[s0.to_a.last+1]
		xs = prev_x + (xs - prev_x)*simplex.solution[0]

		puts "MAXIMIZATION OF T FOR LP WITH SIMPLEX"
		puts "The value of s: #{s}"
		puts "The value of a: #{a}"
		puts "The value of k: #{k}"
		puts "The value of q: #{q}"
		puts "The value of b: #{b}"
		puts "The value of resutl: #{result}"
		puts "The value of xs: #{xs}"
		puts "The value of nus: #{nus}"
		puts "The value of prev_x: #{prev_x}"
		puts "Current status: is feasable: #{check_feasability(a0, xs, b0)}, is optimal: #{check_optimizer(nus)}"
	end
end
