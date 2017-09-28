#lang racket
(require "declarations.rkt")
(provide buildTree)
(provide calcForces)
(provide moveparticles)
(define (inside? area)
  (lambda (p)
    (let* [(x (vec-x (particle-posn p))) 
           (y (vec-y (particle-posn p)))
           (x1 (bbox-llx area))
           (x2 (bbox-rux area))
           (y1 (bbox-lly area))
           (y2 (bbox-ruy area))]
      (cond [(and (> x x1) (< x x2) (> y y1) (< y y2)) #t]
            [(and (= x x1) (>= y y1) (< y y2)) #t]
            [(and (= y y1) (> x x1) (< x x 2)) #t]
            [else #f]))))
(define (present-in-area area particles)
  (filter (inside? area) particles))
(define (avg q w) (/ (+ q w) 2))
(define (LU area)
  (let [(x1 (bbox-llx area))
        (x2 (bbox-rux area))
        (y1 (bbox-lly area))
        (y2 (bbox-ruy area))]
    (bbox x1 (avg y1 y2) (avg x1 x2) y2)))
(define (LB area)
  (let [(x1 (bbox-llx area))
        (x2 (bbox-rux area))
        (y1 (bbox-lly area))
        (y2 (bbox-ruy area))]
    (bbox x1 y1 (avg x1 x2) (avg y1 y2))))
(define (RU area)
  (let [(x1 (bbox-llx area))
        (x2 (bbox-rux area))
        (y1 (bbox-lly area))
        (y2 (bbox-ruy area))]
    (bbox (avg x1 x2) (avg y1 y2) x2 y2)))
(define (RB area)
  (let [(x1 (bbox-llx area))
        (x2 (bbox-rux area))
        (y1 (bbox-lly area))
        (y2 (bbox-ruy area))]
    (bbox (avg x1 x2) y1 x2 (avg y1 y2))))
(define (total-mass particles)
  (define (g a b) (+ (particle-mass a) b))
  (foldr g 0 particles))
(define (centroid particles)
  (define (g a b) (+ (* (particle-mass a) (vec-x (particle-posn a))) b))
  (define (h c d) (+ (* (particle-mass c) (vec-y (particle-posn c))) d))
  (let* ([n (total-mass particles)]
         [X (/ (foldr g 0 particles) n)]
         [Y (/ (foldr h 0 particles) n)])
    (vec X Y)))
(define (buildTree area particles)
  (let* ([lst (present-in-area area particles)]
         [s1 (LU area)]
         [s2 (RU area)]
         [s3 (LB area)]
         [s4 (RB area)]) 
    (cond
      [(> (length lst) 1)
           (gnode (total-mass particles)
                  (centroid particles)
                  (filter gnode? (list (buildTree s1 (present-in-area s1 particles))
                                       (buildTree s2 (present-in-area s2 particles))
                                       (buildTree s3 (present-in-area s3 particles))
                                       (buildTree s4 (present-in-area s4 particles)))))]
      [(= 1 (length lst)) (gnode (particle-mass (car lst))
                                 (particle-posn (car lst))
                                 '())])))
(define (side area)
  (- (bbox-rux area) (bbox-llx area)))
(define (d p C)
  (let* ([P (particle-posn p)]
         [x-p (vec-x P)]
         [y-p (vec-y P)]
         [x-c (vec-x C)]
         [y-c (vec-y C)])
    (sqrt (+ (expt (- x-p x-c) 2) (expt (- y-p y-c) 2)))))

(define (v-add v1 v2)
  (vec (+ (vec-x v1) (vec-x v2)) (+ (vec-y v1) (vec-y v2))))
(define (v-subtract v1 v2)
  (vec (- (vec-x v2) (vec-x v1)) (- (vec-y v2) (vec-y v1))))
(define (far? p s pos)
  (if (> (d p pos) (* theta s)) #t #f))
(define (force p side t)
  (if (equal? (particle-posn p) (gnode-posn t))
      (vec 0 0)
      (let* ([mp (particle-mass p)]
             [mc (gnode-mass t)]
             [location-c (gnode-posn t)]
             [rbar (v-subtract (particle-posn p) location-c)] 
             [r (d p location-c)]
             [common (* g mp mc (/ 1 (expt r 3)))]
             [far (far? p side location-c)])
        (define (f-force lst)
          (if (null? lst) (vec 0 0)
              (v-add (force p (/ side 2) (car lst)) (f-force (cdr lst)))))
        (cond [far
               (vec (* common (vec-x rbar)) 
                    (* common (vec-y rbar)))]
              [(and (not far) (null? (gnode-subtrees t)))
               (vec (* common (vec-x rbar)) 
                    (* common (vec-y rbar)))]
              [else (f-force (gnode-subtrees t))]))))
(define (calcForces area t particles)
  (let ([s (side area)])
  (if (null? particles) '()
      (cons (force (car particles) s t) (calcForces area t (cdr particles))))))
(define (moveparticles particles forces)
  (define (new-state p f)
    (let* ([x (vec-x (particle-posn p))]
          [y (vec-y (particle-posn p))]
          [m (particle-mass p)]
          [vx (vec-x (particle-velocity p))]
          [vy (vec-y (particle-velocity p))]
          [fx (vec-x f)]
          [fy (vec-y f)]
          [ax (/ fx m)]
          [ay (/ fy m)])
      (particle m
                (vec (+ x (+ (* vx timeslice) (* (/ 1 2) ax (expt timeslice 2))))
                     (+ y (* vy timeslice) (* (/ 1 2) ay (expt timeslice 2))))
                (vec (+ vx (* ax timeslice))
                     (+ vy (* ay timeslice))))))
  (if (null? particles) '()
      (cons (new-state (car particles) (car forces))
            (moveparticles (cdr particles) (cdr forces)))))
